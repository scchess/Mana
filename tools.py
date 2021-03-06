import os


def info(x, w):
    x = "[INFO]: " + x
    print(x)
    w.write(x + "\n")


def formatDP(x):
    return '{0:.2f}'.format(x)


def system(cmd, w):
    try:
        info(cmd, w)
        os.system(cmd)
    except Exception:
        pass
