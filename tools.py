import os
import settings


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
