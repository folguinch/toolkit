import datetime

def get_timestamp(fmt=r"%Y%m%d_%H%M%S%f"):
    idf = datetime.datetime.now()
    return idf.strftime(fmt)
