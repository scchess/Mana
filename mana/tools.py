#!/usr/bin/env python
import os
import pickle
import mana.settings as settings


def save(x, file):
    with open(file, "wb") as w:
        pickle.dump(x, w)


def load(file):
    with open(file, "rb") as r:
        return pickle.load(r)


def formatDP(x):
    return '{0:.2f}'.format(x)


def info(x):
    x = "[INFO]: " + x
    print(x)
    with open(settings.LOG_FILE(), "a") as w:
        w.write(x + "\n")


def run(cmd):
    try:
        info(cmd)
        os.system(cmd)
    except Exception:
        pass
