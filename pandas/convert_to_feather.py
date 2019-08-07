from root_pandas import read_root
import sys
from pathlib import Path

filename = sys.argv[1]
feathername = Path(filename).stem + ".feather"
print("converting file ", filename, " to ", feathername)

df = read_root(filename)
df.to_feather(feathername)
