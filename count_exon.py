
#!/usr/bin/env python
#coding=UTF-8

'''
Usage:本脚本用于统计cuffcompare软件输出的gff文件（*.combined.gtf）中对应基因的外显子数目，并根据*.loci文件对其进行分类（PD，GD+PD，GD），
最终绘制不同分类外显子数目的箱线图和折线图()
Author:chenzx(chenzx@biobreeding.com.cn)
Date:2016-10-12
Version:1.0
'''

import argparse
import os
import sys
import matplotlib.pyplot as plt
import numpy as np
import time
global rd
rd = os.getcwd()

def opt():
    args = argparse.ArgumentParser(usage="%(prog)s[options]")
    args.add_argument("--gff",help="gff file for cuffcompare output")
    args.add_argument("--loci",help="loci file for cuffcompare output")
    return args.parse_args()

def ltime():
    return time.strftime("%H:%M:%S",time.localtime(time.time()))

def load_gff(gff):
    '''
    :param gff: cuffcompare输出的gff文件
    :return:一个多级字典
    '''
    sys.stdout.write("[{0}] {1} loading!\n".format(ltime(),os.path.basename(gff)))
    tmpdir = {}
    for line in open(gff,"r"):
        if not line.startswith("#"):
            ele = line.strip().split("\t")
            ids = ele[-1].split(";")
            gene = ids[0].split()[1].strip('"')
            trans = ids[1].split()[1].strip('"')
            exon_num = int(ids[2].split()[1].strip('"'))
            try:
                tmpdir[gene][trans].append(exon_num)
            except KeyError:
                try:
                    tmpdir[gene][trans] = [exon_num]
                except KeyError:
                    tmpdir[gene] = {}
                    tmpdir[gene][trans] = [exon_num]
    var = sys.stdout.write("[{0}] {1} loading successfuly\n".format(ltime(),os.path.basename(gff)))
    return tmpdir

def load_loci(loci):
    '''
    :param loci:cuffcompare输出的loci文件
    :return: 一个嵌套数组的字典
    '''
    sys.stdout.write("[{0}] {1} loading\n".format(ltime(),os.path.basename(loci)))
    tmpdir = {}
    for line in open(loci,"r"):
        if not line.startswith("#"):
            ele = line.strip().split()
            gids = ele[0]
            pd = ele[-2]
            gd = ele[-1]
            if pd == "-" and gd != "-":
                try:
                    tmpdir["GD"].append(gids)
                except KeyError:
                    tmpdir["GD"] = [gids]
            elif pd != "-" and gd == "-":
                try:
                    tmpdir["PD"].append(gids)
                except KeyError:
                    tmpdir["PD"] = [gids]
            elif pd != "-" and gd != "-":
                try:
                    tmpdir["PD+GD"].append(gids)
                except KeyError:
                    tmpdir["PD+GD"] = [gids]
    sys.stdout.write("[{0}] {1} loading successfuly\n".format(ltime(),os.path.basename(loci)))
    return tmpdir

def exon_count(gff,loci):
    '''
    :param gff: cuffcompare output gff file
    :param loci: cuffcompare output loci file
    :return: exon numbers
    '''
    tmpdir = {}
    locidir = load_loci(loci)
    gffdir = load_gff(gff)

    for keys in gffdir.keys():
        if keys in locidir["PD"]:
            for tuple in gffdir[keys].values():
                try:
                    tmpdir["PD"].append(max(tuple))
                except:
                    tmpdir["PD"] = [max(tuple)]
        elif keys in locidir["GD"]:
            for tuple in gffdir[keys].values():
                try:
                    tmpdir["GD"].append(max(tuple))
                except:
                    tmpdir["GD"] = [max(tuple)]
        elif keys in locidir["PD+GD"]:
            for tuple in gffdir[keys].values():
                try:
                    tmpdir["PD+GD"].append(max(tuple))
                except:
                    tmpdir["PD+GD"] = [max(tuple)]
    return tmpdir

def write2file(gff,loci):
    try:
        os.stat("{0}/tmp".format(rd))
    except:
        os.mkdir("{0}/tmp".format(rd))
    tmpdir = exon_count(gff,loci)
    sys.stdout.write("[{0}] Write information to exon_num.txt\n".format(ltime()))
    out = open("{0}/tmp/data4boxplot.txt".format(rd),"w")
    for keys in tmpdir.keys():
        for ids in tmpdir[keys]:
            out.write("{0}\t{1}\n".format(keys,ids))
    sys.stdout.write("[{0}] Write information to exon_num.txt successfuly\n".format(ltime()))

def draw_boxplot():
    try:
        os.stat("{0}/fig".format(rd))
    except:
        os.mkdir("{0}/fig".format(rd))
    try:
        os.stat("{0}/rscript".format(rd))
    except:
        os.mkdir("{0}/rscript".format(rd))

    sys.stdout.write("[{0}] Write R script into boxplot.r\n".format(ltime()))
    out = open("{0}/rscript/boxplot.r".format(rd),"w")
    out.write('library(ggplot2);\n'
              'data <- read.table("{0}/tmp/data4boxplot.txt",header=FALSE,sep="\t");\n'
              'p <- ggplot(data,aes(x=data$V1,y=data$V2,fill=data$V1)) + geom_boxplot() + theme_bw();\n'
              'p + labs(title = "Boxplot for exon numbers") + xlab("Class") + ylab("Exon Numbers") + '
              'theme(axis.title.x = element_text(size = 18,colour = "black"),axis.title.y = element_text(size=18,colour = "black"),'
              'axis.text.x = element_text(size = 18,colour = "black"),axis.text.y = element_text(size=18,colour = "black"),'
              'plot.title = element_text(size=18,color="black"),panel.grid=element_blank(),legend.title=element_text(color="black",size = 18),'
              'legend.text = element_text(color="black",size=18)) + guides(fill = guide_legend(title="Class"));\n'
              'ggsave("{0}/fig/exon_boxplot.png",width=4.5,height = 5);'.format(rd))

    os1 = os.system("Rscript {0}/rscript/boxplot.r".format(rd))
    if os1:
        exit("[{0}] Error in boxplot.r".format(ltime()))
    else:
        sys.stdout.write("[{0}] Draw boxplot successfuly\n")

def draw_line():
    out = open("{0}/tmp/data4line.txt".format(rd),"w")
    tmpdir = {}
    for line in open("{0}/tmp/data4boxplot.txt".format(rd),"r"):
        ele = line.strip().split("\t")
        try:
            tmpdir[ele[0]].append(ele[1])
        except:
            tmpdir[ele[0]] = [ele[1]]
    tmpdir1 = {}
    for keys in tmpdir.keys():
        for ids in sorted(set(tmpdir[keys])):
            try:
                tmpdir1[keys].append(tmpdir[keys].count(ids))
            except:
                tmpdir1[keys] = [tmpdir[keys].count(ids)]

    for keys in tmpdir.keys():
        for ids in sorted(set(tmpdir[keys])):
            per = float(tmpdir[keys].count(ids))/float(sum(tmpdir1[keys]))
            out.write("{0}\t{1}\t{2}\n".format(keys,ids,per))
    sys.stdout.write("[{0}] Write R script to exon_line.r".format(ltime()))
    out1 = open("{0}/rscript/exon_line.r".format(rd),"w")
    out1.write("library(ggplot2);\n"
               "data <- read.table('{0}/tmp/data4line.txt',header = FALSE, sep = '\t');\n"
               "p <- ggplot(data,aes(x = data$V2,y = data$V3, colour = data$V1, group = data$V1)) + geom_line(size = 1);\n"
               "p + theme_bw() + labs(title = 'Line Chart for Exon Numbers') + xlab('Class') + ylab('Percent(%)') + theme(axis.title.x = element_text(size = 18,colour = 'black'),"
               "axis.title.y = element_text(size=18,colour = 'black'),axis.text.x = element_text(size = 18,colour = 'black'),axis.text.y = element_text(size=18,colour = 'black'),"
               "plot.title = element_text(size=18,color='black'),panel.grid=element_blank(),legend.title=element_text(color='black',size = 18),legend.text = element_text(color='black',size=18)) + "
               "guides(color = guide_legend(title='Class'));\n"
               "ggsave('{0}/fig/exon_line.png')".format(rd))
    sys.stdout.write("[{0}] Write R script to exon_line.r successfully\n".format(ltime()))

    sys.stdout.write("[{1}]Rscript {0}/rscript/exon_line.r".format(rd,ltime()))
    os1 = os.system("Rscript {0}/rscript/exon_line.r".format(rd))
    if os1:
        exit("[{0}] Error in exon_line.r\n".format(ltime()))
    else:
        sys.stdout.write("[{0}] Draw line successfully\n".format(ltime()))
    sys.stdout.write("[{0}] Draw line successfuly\n".format(rd))

def main():
    ags = opt()
    if ags.gff == None or ags.loci == None:
        os.system("python {0} -h".format(sys.argv[0]))
        exit("[{0}] Error Incomplete Options\n".format(ltime()))
    gff = ags.gff
    loci = ags.loci
    write2file(gff,loci)
    draw_boxplot()
    draw_line()

if __name__ == "__main__":
    main()
