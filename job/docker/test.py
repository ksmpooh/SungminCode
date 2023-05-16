############## CODE 1
#### test file IO

df = open("filein.txt","r")
out = open("fileout.txt","w")

while 1:
    line = df.realine() # line 읽기
    if not line: # line 없으면 break
        print("Not line.. done")
        break
    # 필요한 잡업 수행

    #
    out.write("%s\n") # line 작업 완료되자마자 바로 write

out.close()

############## CODE 2 : test 용
#### line 4개만 할 때 테스트 , count 추가

df = open("filein.txt","r")
out = open("fileout.txt","w")
count = 0
while 1:
    count = count + 1
    line = df.realine() # line 읽기
    if count == 10: # count 10이면 break! 
        print("Count : %s.. done"%str(count))
        break
    if not line: # line 없으면 break
        print("Not line.. done")
        break
    # 필요한 잡업 수행

    #
    out.write("%s\n") # line 작업 완료되자마자 바로 write

out.close()