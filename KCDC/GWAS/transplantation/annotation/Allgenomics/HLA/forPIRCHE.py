import sys

# 입력 파일과 출력 파일 이름 설정
input_file = sys.argv[1]
output_file = input_file.replace(".csv", ".processing.csv")

# 파일 읽고 쓰기 (with 사용하여 자동으로 닫힘)
with open(input_file, "r") as df, open(output_file, "w") as out:
    out.write(df.readline())  # 첫 번째 줄 복사

    while True:
        line = df.readline()
        if not line:  # 파일 끝이면 종료
            break
        out.write(line)

        line = df.readline()
        if not line:  # 파일 끝이면 종료
            break
        out.write(line)

        out.write("\n")  # 두 줄마다 공백 줄 추가