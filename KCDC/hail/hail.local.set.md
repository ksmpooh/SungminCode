# (Local) hail 설치 (in conda env) : 20240103

## 1. 설치

### 사전준비 (필요한 package 설치)

apt-get install openjdk-11-jre-headless #JAVA 11로 셋팅
apt-get install libopenblas-base liblapack3 #추가 필요 library

### conda 가상환경 만들기
conda create -n [가상환경이름] python=3.9 #최신 hail에서 3.9이상 버전 권장

conda create -n hail python=3.9 #hail 이라는 이름을 가진 python3.9 버전의 conda 가상환경 만들기

### conda 가상환경 실행 및 hail 설치
conda env list #conda 환경 확인
conda activate hail #hail 이라는 가상환경 실행
python -m pip install hail #python에 hail 설치

## 2. hail 실행
### conda 가상환경 확인
conda env list

### hail 가상환경 시작
conda activate hail

### python 실행, import hail
python
import hail as hl
