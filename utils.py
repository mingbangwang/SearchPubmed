#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import pandas as pd


def dataframe_dir_to_df(dataframe_dir):
    """
    导入dataframe_dir，里面每一个df的表头是一样的，将之合并为一个dataframe，并返回列表
    """
    df_files = os.listdir(dataframe_dir)
    # print(df_files)
    df = pd.DataFrame()
    for df_file in df_files:
        df_pubmed = pd.read_csv(os.path.join(dataframe_dir, str(df_file)), sep="\t", index_col=False)
        df = pd.concat([df, df_pubmed], axis=0)
    return df