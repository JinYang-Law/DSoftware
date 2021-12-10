# -*- coding: UTF-8 -*-
import logging
import numpy as np
import os
import argparse
import sys
import csv
import pandas as pd
from collections import defaultdict
from itertools import groupby
from operator import itemgetter
import time
from functools import wraps


def timefn(fn):
    """计算性能的修饰器"""

    @wraps(fn)
    def measure_time(*args, **kwargs):
        t1 = time.time()
        result = fn(*args, **kwargs)
        t2 = time.time()
        print(f"@timefn: {fn.__name__} took {t2 - t1: .5f} s")
        return result

    return measure_time


class TmbCalculation(object):
    """
        利用暴力枚举的方法进行操作
    """

    def __init__(self):
        pass

    def _ValueSorted(self, inFile):
        """
        :param inFile: bed 文件， 要求已经 按照 chr - start - end 排序 sort -k1,1V -k2,2n -k3,3n bed.txt
        :return:
        """
        region = 0
        bedDict = defaultdict(list)

        with open(inFile, newline='') as csvfile:
            reader = csv.DictReader(csvfile, fieldnames=['Chr', 'Start', 'End'], delimiter="\t")
            for row in reader:
                bedDict[row['Chr']].append((int(row['Start']), int(row['End'])))
                region += int(row['End']) - int(row['Start']) + 1
                # bedDict = sorted(bedDict.items(), key=lambda x:x[1])
        for k, v in bedDict.items():
            bedDict[k] = sorted(v, key=lambda x: x[0])

        return region, bedDict

    def _sorted(self, inFile):
        """
        :param inFile: 黑名单： vcf 格式
        :return: 排序的数组
        """
        bList = []
        with open(inFile, 'r') as csvfile:
            spamreader = csv.reader(csvfile, delimiter='\t', quotechar='|')
            for row in spamreader:
                bList.append('_'.join(row[:5]))
        bList = sorted(bList)

        return bList

    @timefn
    def tmbCalc(self, varfile, outfile, bedfile, minfreq, cutoff, blacklist):
        # 构建 bed 字典， chrom : {(start, end)...} ;并按value 排序
        region, bedDict = self._ValueSorted(bedfile)

        # 构建黑名单数组，并排序， 构建字典， 直接get() 会不会更快？
        bList = []
        if blacklist:
            bList = self._sorted(blacklist)

        # 处理变异文件
        totalNum = 0
        with open(varfile, newline='') as csvfile:
            reader = csv.DictReader(csvfile, delimiter="\t")
            for row in reader:
                if not bedDict.get(row['Chr']):
                    continue
                elif float(row['Freq']) < minfreq:
                    continue
                elif '_'.join(list(map(lambda x: row[x], ["Chr", "Start", "End", "Ref", "Alt"]))) in bList:
                    continue

                for pos in bedDict[row['Chr']]:
                    if int(row['Start']) > pos[1]:
                        continue
                    elif int(row['End']) < pos[0]:
                        break
                    else:
                        totalNum += 1

        tmb = totalNum / region * 1000000

        if tmb >= cutoff:
            tmbStatus = 'TMB-H'
        else:
            tmbStatus = 'TMB-L'

        print("Target Region : %s\n" % region)
        print("Total Mutation : %s\n" % totalNum)
        print("TMB Value : %s\n" % tmb)
        print("TMB Status : %s\n" % tmbStatus)

        with open(outfile, "w") as wf:
            wf.writelines("Target Region : %s\n" % region)
            wf.writelines("Total Mutation : %s\n" % totalNum)
            wf.writelines("TMB Value : %s\n" % tmb)
            wf.writelines("TMB Status : %s\n" % tmbStatus)
        pass


class TmbSolution(object):
    """
        利用 贪心算法进行操作
    """

    def __init__(self):
        pass

    def _sorted(self, inFile):
        """
        :param inFile: 黑名单： vcf 格式
        :return: 排序的数组
        """
        bList = []
        with open(inFile, 'r') as csvfile:
            spamreader = csv.reader(csvfile, delimiter='\t', quotechar='|')
            for row in spamreader:
                bList.append('_'.join(row[:5]))
        bList = sorted(bList)

        return bList


    def _ValueSorted(self, inFile):
        """
        :param inFile: bed 文件， 要求已经 按照 chr - start - end 排序 sort -k1,1V -k2,2n -k3,3n bed.txt
        :return:
        """
        region = 0
        bedDict = defaultdict(list)

        with open(inFile, newline='') as csvfile:
            reader = csv.DictReader(csvfile, fieldnames=['Chr', 'Start', 'End'], delimiter="\t")
            for row in reader:
                bedDict[row['Chr']].append((int(row['Start']), int(row['End'])))
                region += int(row['End']) - int(row['Start']) + 1
                # bedDict = sorted(bedDict.items(), key=lambda x:x[1])
        for k, v in bedDict.items():
            bedDict[k] = sorted(v, key=lambda x: x[0])

        return region, bedDict

    @timefn
    def tmbCalc(self, varfile, outfile, bedfile, minfreq, cutoff, blacklist):

        # 构建黑名单数组，并排序
        bList = []
        if blacklist:
            bList = self._sorted(blacklist)

        # 构建 bed 字典， chrom : {(start, end)...} ;并按value 排序
        region, bedDict = self._ValueSorted(bedfile)

        # 处理变异文件
        with open(varfile, newline='') as csvfile:
            reader = csv.DictReader(csvfile, delimiter="\t")
            info = sorted(reader, key=lambda x: (x.get("Chr"), int(x.get("Start")), int(x.get("End"))))
            totalNum = 0

            for ind, group in groupby(info, key=itemgetter('Chr')):

                if not bedDict.get(ind):
                    continue
                group = list(group)
                varCur = 0
                bedCur = 0

                while varCur < len(group) and bedCur < len(bedDict[ind]):

                    if float(group[varCur]['Freq']) < minfreq or '_'.join(
                            list(map(lambda x: group[varCur][x], ["Chr", "Start", "End", "Ref",
                                                                  "Alt"]))) in bList:  # 这里是O(n) 其实改，1. 大数据构建字典-慢 2. in O(n) 3.数据不多没关系
                        varCur += 1
                        continue

                    if int(group[varCur]['Start']) > int(bedDict[ind][bedCur][1]):
                        bedCur += 1
                        continue
                    elif int(group[varCur]['End']) < int(bedDict[ind][bedCur][0]):
                        varCur += 1
                        continue
                    else:
                        totalNum += 1
                        varCur += 1

        tmb = totalNum / region * 1000000

        if tmb >= cutoff:
            tmbStatus = 'TMB-H'
        else:
            tmbStatus = 'TMB-L'

        print("Target Region : %s\n" % region)
        print("Total Mutation : %s\n" % totalNum)
        print("TMB Value : %s\n" % tmb)
        print("TMB Status : %s\n" % tmbStatus)

        with open(outfile, "w") as wf:
            wf.writelines("Target Region : %s\n" % region)
            wf.writelines("Total Mutation : %s\n" % totalNum)
            wf.writelines("TMB Value : %s\n" % tmb)
            wf.writelines("TMB Status : %s\n" % tmbStatus)
        pass


class TmbCalcPandas(object):
    """
     利用 pandas 进行操作
    """

    def __init__(self):
        pass

    def bed2interval(self, bed):
        """
        bed : 1-base
        :param bed: bed 文件
        :return: region && interval pandas 对象
        """
        chrom = []
        pos = []
        region = 0
        with open(bed, newline='') as csvfile:
            spamreader = csv.reader(csvfile, delimiter='\t')
            for row in spamreader:
                region += int(row[2]) - int(row[1]) + 1
                pos.extend([i for i in range(int(row[1]), int(row[2]) + 1)])
                chrom.extend([row[0]] * (int(row[2]) - int(row[1]) + 1))

        bedDataF = pd.DataFrame({'Chr': chrom, 'Start': pos, 'End': pos})
        return region, bedDataF

    def read_file(self, input):
        """
        :param input: 变异检测文件， 列至少包括：Chr	Start	End	Ref	Alt Freq
        :return: pandas 对象
        """
        return pd.read_csv(input, sep="\t", dtype={'Freq': np.float64})

    @timefn
    def tmbCalc(self, varfile, outfile, bedfile, minfreq, cutoff, blacklist):
        # 读取 变异检测文件
        varData = self.read_file(varfile)

        # 读取 panel 文件
        region, bedData = self.bed2interval(bedfile)

        # 筛选符合条件的突变
        onBed = pd.merge(varData, bedData, how='inner', on=["Chr", 'Start'], suffixes=("_var", "_bed"))
        gtFreqBed = onBed[(onBed['Freq'] >= minfreq)]

        # 过滤黑名单的突变
        if blacklist:
            blistData = pd.read_csv(blacklist, sep="\t")
            gtFreqBed = pd.merge(gtFreqBed, blistData, how='left', on=["Chr", 'Start', 'Ref', 'Alt'],
                                 indicator=True, suffixes=("_var", "_filter"))
            gtFreqBed = gtFreqBed[gtFreqBed['_merge'] == "left_only"]

        tmb = gtFreqBed.shape[0] / region * 1000000

        if tmb >= cutoff:
            tmbStatus = 'TMB-H'
        else:
            tmbStatus = 'TMB-L'

        with open(outfile, "w") as wf:
            wf.writelines("Target Region : %s\n" % region)
            wf.writelines("Total Mutation : %s\n" % gtFreqBed.shape[0])
            wf.writelines("TMB Value : %s\n" % tmb)
            wf.writelines("TMB Status : %s\n" % tmbStatus)
        pass


def getParser():
    # parser = argparse.ArgumentParser(description="")   ## 创建一个解析对象，description = 描述内容
    # parser.add_argument('-v',"--verbosity",type=int,choices=[0,1,2],\
    #                 default=1,required=True,help="increase output verbosity")   ## 向对象添加命令行参数和选项
    # parser.parse_args()     ## 解析

    # required
    parser = argparse.ArgumentParser(description="calculate TMB")
    parser.add_argument('-i', '--var', dest="input", type=str, required=True, help="")
    parser.add_argument('-o', '--out', dest="output", type=str, required=True, help="")
    parser.add_argument('-b', '--bed', dest="bedfile", type=str, required=True, help="")

    # optional
    parser.add_argument('-f', '--freq', dest="minfreq", type=float, required=False, default=0.05,
                        help="[default: %(default)s]")
    parser.add_argument('-bl', '--blist', dest="blacklist", type=float, required=False, help="[optional: %(default)s]")
    parser.add_argument('-cf', '--coff', dest="cutoff", type=int, required=False, default=10,
                        help="[default: %(default)s]")

    # 返回对象
    return parser.parse_args()


def main():
    args = getParser()
    input = args.input
    output = args.output
    bedfile = args.bedfile
    minfreq = args.minfreq
    blacklist = args.blacklist
    cutoff = args.cutoff


    # 对输入参数进行判断
    for file in [input, bedfile]:
        if not os.path.exists(file):
            logging.error("ERROR: file named %s does not exists !!!" % file)
            sys.exit()

    # 实例化 TMB 检测 对象
    # tmbObject = TmbCalculation()
    # tmbObject.tmbCalc(input, output, bedfile, minfreq, cutoff, blacklist)

    tmbObject = TmbSolution()
    tmbObject.tmbCalc(input, output, bedfile, minfreq, cutoff, blacklist)

    # tmbObject = TmbCalcPandas()
    # tmbObject.tmbCalc(input, output, bedfile, minfreq, cutoff, blacklist)


if __name__ == '__main__':
    main()
