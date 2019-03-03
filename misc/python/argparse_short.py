import argparse

parser = argparse.ArgumentParser(description='rebuild matrix options')

parser.add_argument('-v', action="store_true", default=False)
parser.add_argument('-csv', action="store", dest="csv")
parser.add_argument('-h5', action="store_true", default=False)
#parser.add_argument('-npx', action="store", dest="c", type=int)

#print(parser.parse_args(['-a', '-bval', '-c', '3']))
args = parser.parse_args()

print(args.v)
print(args.csv)
print(args.h5)

