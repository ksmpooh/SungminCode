### 두번째 시간

y <- c(1,2,3,NA)
y
is.na(y)
summary(y)

ages <- c(48,78,56,88,-999,13,26,-999)
ages
ages[ages == -999] = NA
ages

sum(is.na(ages))
is.na(ages)

is.na(NA)
is.na("NA")

sum(ages) # 결측치가 있으면 합계가 안구해진다. 평균도 안된다.
 
#결측치를 있을경우 옵션 추가하여 구한다.

sum(ages,na.rm = TRUE)
mean(ages,na.rm  = TRUE)
                                              
weight <- c(65.4,55,380,72.2,51,NA)       
height <-c(170,155,NA,173, 161, 166)                                                               
gender <- c("M","F","M","M","F","F")

#gender <- c(1,2,1,1,2,2)
#gender <- factor(gender,levels = c(1,2),labels = c("M","F"))  이거 좋은듯
#gender

testDate <- c("2013/09/01","2013/09/01","2013/09/05","2013/09/14","2013/10/11","2013/10/26")

patients <- data.frame(weight = weight, height = height, gender = gender, testDate = testDate)
patients

str(patients)
na.omit(patients)
complete.cases(patients)
patients[complete.cases(patients),]

patients.sub <- patients[,c("weight","height")]
patients.sub

apply(patients.sub,2,mean,na.rm = TRUE) # 1 = row, 2 = col  ##함수를 지속적으로 계산할때 apply를 사용하여 계산다.


#날짜를 계산할수있게 만든다.
as.Date(testDate)
patients
class(as.Date(testDate))

patients$testDate <- as.Date(testDate)

patients$testDate[5] - patients$testDate[1]
today<- Sys.Date()
today

difftime(today,patients$testDate[1],units = "weeks") # 주단위로.. 선택하여 계산
as.numeric(difftime(today,patients$testDate[1],units = "weeks"))

a <- c(4,6,8)
a
is.numeric(a)
is.matrix(a)
a
a <- as.character(a)
a
a <- as.factor(a)
a

for(i in 1:5){
  print(i)
  print("Hello R")
}

i <- 5
while( i > 0){
  print("Hello")
  i <- i - 1
}

#조건문
ifelse(patients$weight > 150, NA, patients$weight)
patients$weight[patients$weight>150] 
patients$weight <-ifelse(patients$weight > 150, NA, patients$weight) #변수에 저장하여 계산한다.
patients

scores <- c(55,70,80,90,59,87)
scores
ifelse(scores>= 60, "passed","Fail")

isHeightOver160 <-ifelse(patients$height > 160, "over 160","under 160")
isHeightOver160
patients <- cbind(patients,data.frame(isHeightOver160 = isHeightOver160))
patients

lapply(patients.sub,mean,na.rm = TRUE)
sapply(patients.sub,mean,na.rm = TRUE)

#tapply 각각 범주에 대해서.. apply
tapply(patients$weight,patients$gender, mean, na.rm = TRUE)
tapply(patients$height,patients$gender, mean, na.rm = TRUE)


#subset
patients[patients$weight>60, patients$gender == "M"]
options(width = 240)
patients[(patients$weight>60 & patients$gender == "M" & !is.na(patients$weight) & !is.na(patients$gender)),]
subset(patients,weight > 60 & gender == "M")
subset(patients,weight > 60 & gender == "M",select = c("weight","height"))


## transform
cbind(patients,data.frame(weight_pounds = 2.2* patients$weight))
transform(patients,weight_pounds = 2.2*weight)
transform(patients,weight_pounds = 2.2*weight,height_inch = 0.39*height)

#aggregate
aggregate(patients[,1:2],by = list(gender, isHeightOver160),mean,na.rm = TRUE)

#table
table(patients$gender,patients$isHeightOver160)
