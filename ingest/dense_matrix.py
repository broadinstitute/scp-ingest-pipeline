

class Dense():
    def __init__(self, file_path):
        if not os.path.exists(file_path):
		          raise IOError(f"File '{file_path}' not found")
        self.file = open(file_path,'r')
        self.file_name, self.filetype = os.path.splitext(file_path)

def extract(file, size=500):
    #skip first line for dense matrix
    self.file.readline()
    while True:
        next_lines = list(islice(file, number_of_lines))
        if not next_lines:
            # Let's worker function that there are no more extracted data
            return []
            break
        yield next_lines
