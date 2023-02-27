import sys

def main():
    if not sys.argv[1]:
        a = "defualt"
    else:
        a = sys.argv[1]
    print(a)

main()