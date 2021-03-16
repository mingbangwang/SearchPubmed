#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from Bio import Entrez
from Bio import Medline
from wordcloud import WordCloud
import os
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

# 导入自定义的类或者函数
from utils import *

import warnings

warnings.filterwarnings('ignore')

# 设置工作文件夹
full_path = os.path.realpath(__file__)
work_dir = os.path.dirname(os.path.realpath(__file__))
print(f"Change CWD to: {work_dir}")

# 设置输出文件夹

## 设置项目结果保存目录
out_dir = os.path.join(work_dir, "pubmed_out")
os.makedirs(out_dir, mode=0o777, exist_ok=True)
print("out_dir:", out_dir)

## 设置待查杂志
# 定义一个杂志和ISSN的字典

issn_journal_dict = {
    # nature系列
    "1476-4687": "nature", "2041-1723": "nature_communication", "1097-6256": "nature_neuroscience",

    # science系列
    "0036-8075": "science", "2375-2548": "science_advance", "1946-6234": "science_translational_medicine",

    # cell系列
    "0092-8674": "cell", "2211-1247": "cell_report",

    # nejm
    "0028-4793": "nejm", "1546-170X": "nature_medicine",

    # lancet系列
    "1474-547X": "lancet", "1474-4422": "Lancet_neurology",

    # bmj
    "1756-1833": "bmj",

    # jama系列
    "0098-7484": "jama", "2168-622X": "jama_psychiatry",

    # 神经精神科学
    "0896-6273": "neuron", "1359-4184": "molecular_psychiatry", "0006-3223": "biological_psychiatry", "0006-8950": "brain", "0028-3878": "neurology",

    # 遗传
    "1061-4036": "nature_genetics", "0002-9297": "ajhg", "1530-0366": "genetics_in_medicine",

    # 儿科
    "2168-6203": "jama_pediatrics", "0031-4005": "pediatrics",

    # 医学综合
    "0021-9738": "j_clin_invest", "1549-1277": "plos_medicine",

    # 微生物组相关
    "0017-5749": "gut", "0016-5085": "gastroenterology", "2049-2618": "microbiome", "1931-3128": "cell_host_microbe", "1949-0976": "gut_microbes",

    # 免疫传染病
    "1097-4180": "immunity",

    # 组学杂志
    "1088-9051": "genome_research", "1474-760X": "genome_biology", "1756-994X": "genome_medicine", "2050-084X": "elife", "1672-0229": "genomics_proteomics_bioinformatics",

    # 其他综合
    "2198-3844": "advanced_science", "1097-2765 ": "molecular_cell", "0027-8424": "pnas", "1001-0602": "cell_research", "2405-4712": "cell_systems"
}

# # 测试获取杂志的ISSN
# issn_journal_dict.get("2041-1723")
# issn_journal_dict['2041-1723']


# 获取所有杂志的issn，并形成一个字符串以“or"链接
issn_terms = []
journal_list = []
for key in issn_journal_dict.keys():
    #     print(key)
    issn_terms.append(key)
    journal_list.append(issn_journal_dict[key])
# print(issn_terms)
print("journal_list:", journal_list)

## 设置搜索term
keyword_terms = ["deep learning", "machine learning", "neural network"]
print("keyword_terms:", keyword_terms)

# keyword_term = keyword_terms[0]  # 可以设置循环
for keyword_term in keyword_terms:
    print("pubmed for keyword_term:", keyword_term)
    keyword_term_save = keyword_term.replace(" ", "_")

    # 建立保存原始keyword_term pubmed原始结果的文件夹
    keyword_term_dir = os.path.join(out_dir, keyword_term_save)
    keyword_term_rawdata = os.path.join(keyword_term_dir, "pubmed_tsv")

    os.makedirs(keyword_term_rawdata, mode=0o777, exist_ok=True)

    ## 开始进行pubmed查找
    # 针对指定的keyword_term查询每一个杂志
    # issn_term = issn_terms[1] # 可以设置循环
    for issn_term in issn_terms:
    # print(issn)
        print("searching {} in {}".format(keyword_term,issn_journal_dict[issn_term]))
        # 设置keyword term
        # term = "deep learning and {}".format(issn) # 2041-1723,nature communication等
        term = keyword_term + " and " + issn_term
        term_list = term.split(" ")
        term_issn_save = keyword_term_save + "_" + issn_journal_dict[issn_term]
        # term_issn_save

        # search pubmed
        # 如果已经跑过这个杂志，就跳过
        out_file = os.path.join(keyword_term_rawdata, "pubmed_out_{}.tsv".format(term_issn_save))
        if os.path.exists(out_file):
            print("already done: searching for {} in {}".format(keyword_term,issn_journal_dict[issn_term]))
            continue

        # 开始跑pubmed
        Entrez.email = 'mingbang.wang.bgi@qq.com'  # always tell who you are
        handle = Entrez.egquery(term=term)

        record = Entrez.read(handle)
        for row in record["eGQueryResult"]:
            if row["DbName"] == "pubmed":
                print("count:", row["Count"])

        handle = Entrez.esearch(db="pubmed", term=term, retmax=500)
        record = Entrez.read(handle)
        idlist = record["IdList"]
        # print (idlist[0:40]) #取前40个看下输出的信息是什么

        # id = idlist[0:10] #取前10个测试一下文章标题和作者信息
        id = idlist
        print("records counts:", len(id))

        # 如果没有查到结果，跳到下一个循环
        if len(id) == 0:
            print("none count: searching for {} in {}".format(keyword_term,issn_journal_dict[issn_term]))
            continue

        fout = open(out_file, 'w')
        fout.write("Authors\tTitle\tJournal\tYear\tPublication_citation\tAbstract\tKeywords\tDate_epublished\tPMID\n")

        handle = Entrez.efetch(db="pubmed", id=id, rettype="medline", retmode="text")
        records = Medline.parse(handle)
        records = list(records)  # records 是一个迭代器，所以只能访问这些records一次。如果想保存这些records，需要把他们转成列表。
        journal_list = []
        abstract_list = []
        keyword_list = []
        for record in records:
            Authors_full_name = "; ".join(record.get("FAU", "?"))
            Authors_simplified_name = "; ".join(record.get("AU", "?"))
            Publication_title = record.get("TI", "?")
            Publication_citation = record.get("SO", "?")
            Journal_Title_simplified = record.get("TA", "?")
            Journal_Title = record.get("JT", "?")
            Abstract = record.get("AB", "?")

            # for Keywords
            Keywords_ = record.get("OT", "?")
            # remove "*" in Keywords_
            Keywords = []
            for word in Keywords_:
                Keywords.append(word.replace("*", ""))
            if Keywords == ['?']:
                Keywords = ""
            else:
                Keywords = "; ".join(Keywords)

            PMID = "https://pubmed.ncbi.nlm.nih.gov/" + record.get("PMID", "?")
            Date_epublished = record.get("DP", "?")
            Year_epublished = Date_epublished.split(" ")[0]
            Date_epublished_num = record.get("DEP", "?")
            ISSN = record.get("IS", "?")
            Library_ID = record.get("LID", "?")

            # for word_cloud
            journal_list.append(Journal_Title_simplified)
            abstract_list.append(Abstract)
            keyword_list.append(Keywords)

            fout.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(Authors_full_name, Publication_title, Journal_Title_simplified, Year_epublished, Publication_citation, Abstract, Keywords, Date_epublished_num, PMID))
        fout.close()


    ## 结果统计
    # 调用函数原始的df
    df = dataframe_dir_to_df(keyword_term_rawdata)
    # 只保留research paper
    df = df[df['Abstract'] != '?']  # 可以修改为其他条件
    # 保存df
    df.to_csv(os.path.join(keyword_term_dir, "{}_summary.tsv".format(keyword_term_save)), sep="\t", index=False)

    # 按照杂志进行统计
    result = pd.value_counts(df['Journal'])
    # type(result)
    df_summary = pd.DataFrame(columns=['Journal', 'Count'])
    df_summary['Count'] = result
    df_summary["Journal"] = df_summary.index
    # df_summary
    df_summary = df_summary[df_summary["Count"] > 1]  #
    # df_summary
    plt.figure(figsize=(10, 12))  #
    ax = sns.barplot(y=df_summary.Journal, x=df_summary.Count)  #
    plt.xticks(fontsize=12)
    plt.yticks(fontsize=10)  #
    plt.title('Pubmed term:\n{}'.format(keyword_term), fontsize=14)  # keyword_term 替换project_id
    plt.xlabel("Pubmed records count", fontsize=14)
    plt.ylabel("Journal", fontsize=14)

    plt.show()
    fig = ax.get_figure()
    fig.savefig(keyword_term_dir + "/barplot_{}_journal.pdf".format(keyword_term_save), dpi=400, bbox_inches='tight')
    # df_summary.to_csv("{}_{}_summary.tsv".format(keyword_term,'Journal'),sep="\t",index=False)

    # 按照年份进行统计
    result = pd.value_counts(df['Year'])
    # type(result)
    df_summary = pd.DataFrame(columns=['Year', 'Count'])
    df_summary['Count'] = result
    df_summary["Year"] = df_summary.index
    # df_summary
    df_summary = df_summary[df_summary["Count"] > 1]  #
    # df_summary
    plt.figure(figsize=(10, 12))  #
    ax = sns.barplot(x=df_summary.Year, y=df_summary.Count)  # ,orient='h'
    plt.xticks(fontsize=12)
    plt.yticks(fontsize=10)  #
    plt.title('Pubmed term:\n{}'.format(keyword_term), fontsize=14)  # keyword_term 替换project_id
    plt.ylabel("Pubmed records count", fontsize=14)
    plt.xlabel("Year", fontsize=14)
    # 展示图片
    plt.show()
    # 保存图片
    fig = ax.get_figure()
    fig.savefig(keyword_term_dir + "/barplot_{}_year.pdf".format(keyword_term_save), dpi=400, bbox_inches='tight')
    # df_summary.to_csv("{}_{}_summary.tsv".format(keyword_term,'Journal'),sep="\t",index=False)

    ## 词云图
    df = pd.read_csv(os.path.join(keyword_term_dir, "{}_summary.tsv".format(keyword_term_save)), sep="\t")
    wordcloud_list = ["Title", "Abstract", "Keywords"]  # 有时候提取不出来
    for name in wordcloud_list:
        # remove empty text
        text_list = df[name].to_list()
        text_list_filtered = []
        for text in text_list:
            if str(text) != 'nan':
                text_list_filtered.append(str(text))
        # print(text_list_filtered)
        # word cloud
        text = ' '.join(text_list_filtered)
        wc = WordCloud(width=1200, height=800)
        wc.generate(text)
        wc.to_file(keyword_term_dir + "/wordcloud_{}_{}.pdf".format(keyword_term_save, name))
