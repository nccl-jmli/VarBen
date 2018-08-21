import logging


class InvalidLog(object):
    def __init__(self, log_file):
        self.log_file = log_file
        self.log = open(log_file, 'w')

    def info(self, *args):
        print args
        args_str = ", ".join(args) + "\n"
        self.log.write(args_str)

    def close(self):
        self.log.close()


class RunLog(object):
    def __init__(self, log_file):

        log_format = logging.Formatter(
            "[%(asctime)s] - %(levelname)s - %(name)s - %(pathname)s line:%(lineno)d - %(message)s",
            "%Y-%m-%d %H:%M:%S")
        log_file_handler = logging.FileHandler(log_file)
        log_file_handler.setFormatter(log_format)
        self.log = logging.getLogger()
        self.log.addHandler(log_file_handler)
        self.log.setLevel(logging.DEBUG)

    def info(self, *args):
        args_str = ", ".join(args)
        self.log.info(args_str)


def trace_func(func):
    """
    A decorate function to track all function invoke information with DEBUG level
    Usage:
    @trace_func
    def any_function(any parameter)
    """
    def trace_log(*args, **kargs):
        print 'Start %s(%s, %s)...' % (func.__name__, args, kargs)
        return func(*args, **kargs)
    return trace_log
