#!/bin/bash
# f1: cel_files
# f2: rm cel ID
# 명령어 인수로 파일 경로를 받음
f1=$1
f2=$2

# f2에 있는 패턴들을 배열로 읽어들임
mapfile -t patterns < "$f2"

# f1에서 f2의 패턴을 모두 제거한 결과를 저장할 새로운 파일 생성
output="2nd.cel_files.txt"

# f1 파일에서 각 줄을 확인하며 f2의 패턴이 포함되지 않은 줄만 출력
while IFS= read -r line; do
    # f2에 있는 패턴을 모두 확인하여 포함되지 않으면 출력
    skip=false
    for pattern in "${patterns[@]}"; do
        if [[ "$line" == *"$pattern"* ]]; then
            skip=true
            break
        fi
    done
    # skip이 false일 때만 라인을 출력 (즉, 패턴이 포함되지 않은 경우)
    if [ "$skip" = false ]; then
        echo "$line" >> "$output"
    fi
done < "$f1"

echo "Filtered file created: $output"