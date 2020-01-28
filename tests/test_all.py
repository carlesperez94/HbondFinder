import pytest
import os
import shutil
import pandas as pd
import subprocess

OUT_REP = "hbond_analysis/report_hb_1"
CURR_DIR = os.path.abspath(os.curdir)


def clean():
    if os.path.exists("hbond_analysis"):
        shutil.rmtree("hbond_analysis")


def test_1():
    clean()
    subprocess.call("bash {}/simple_test.sh".format(CURR_DIR).split())
    assert os.path.exists(OUT_REP)


def test_2():
    report = pd.read_csv(OUT_REP, sep="\t", header=0)
    hbond_count = report["92O-94O_H"][0]
    assert hbond_count == 2.0


def test_3():
    clean()
    subprocess.call("bash {}/long_close_test.sh".format(CURR_DIR).split())
    assert os.path.exists(OUT_REP)


def test_4():
    report = pd.read_csv(OUT_REP, sep="\t", header=0)
    hbond_count = float(report["92O-94O_H"][0])
    assert hbond_count == 3.0
