

def copy_fields(destination, source, fields):
    for field in fields:
        destination[field] = source[field]