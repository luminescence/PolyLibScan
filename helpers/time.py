from datetime import datetime
from datetime import timedelta

def date_string():
	return datetime.strftime(datetime.now(), format='%Y-%m-%d')

def time_string():
	return datetime.strftime(datetime.now(), format='%Y-%m-%dT%H:%M:%S')

def add_to_time(time, seconds):
    return time + timedelta(seconds=seconds)

def conversion(time_obj):
    if isinstance(time_obj, basestring):
        if len(time_obj) == 10:
            return datetime.strptime(time_obj, '%Y-%m-%d')
        elif len(time_obj) == 19:
            return datetime.strptime(time_obj, '%Y-%m-%dT%H:%M:%S')
        else:
            raise Exception('Wrong Format.')
    elif isinstance(time_obj, datetime):
        return datetime.strftime(time_obj, format='%Y-%m-%dT%H:%M:%S')
    else:
        raise Exception('Wrong Type.')