def write_txt(text, filename, mode='w'):
    with open(filename, mode) as out:
        out.write(text)
