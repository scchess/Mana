import os
import pickle


def info(x, w):
    x = "[INFO]: " + x
    print(x)
    w.write(x + "\n")


def formatDP(x):
    return '{0:.2f}'.format(x)


def system(cmd, w):
    info(cmd, w)
    os.system(cmd)


def save(x, file):
    with open(file, "wb") as w:
        pickle.dump(x, w, protocol=4)


def load(file):
    with open(file, "rb") as r:
        return pickle.load(r)
