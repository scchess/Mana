import os


def formatDP(x):
    return '{0:.2f}'.format(x)


def system(cmd, logger):
    msg = "[INFO]: " + cmd
    if logger:
        logger.write(msg + "\n")
    print(msg)
    os.system(cmd)
