#3 처음

m <- matrix(1:20, nrow = 5, ncol = 4)
m

m <- matrix(1:20,5,4,byrow = TRUE)
m[1,1]
m[1,]
m[c(3,4),]
m[4,2] <- 0
m

m[5,] <- c(117,118,119,120)
m

m[5,] <- c(1,2) # multiple이 맞으면 된다.
m
m[5,] <- c(1,2,3) # error

m <- matrix(1:20,nrow = 5)
m
ls() # 모든 변수를 확인할때

sum(m)
sum(m[2,])

# apply 는 연속해서 연산을(함수)를 사용할 떄 사용한다.
#apply(x,margin,fun) margin = 1 은 row, 2는c column
apply(m,1,sum)
apply(m,2,sum)

A <- matrix(1:2,2)
B <- matrix(c(4,10),2)
A
B
C <- matrix(c(4,10),1)
C
A %*% C  # matrix 곱셈

m <- matrix(1:20,nrow = 5)
m[m>10]
cbind(m,c(1:5))
m <- cbind(m,c(1:5))
m <- rbind(m,c(10:14))
m

vector <- c(1:5)
vector
names(vector) <- c("A","B","C","D","E")
vector
vector[2]
vector["B"]
names(vector)


countTable <- matrix(c(189,10845,104,10933),nrow =2, byrow =T)
countTable 
rownames(countTable) <-  c("Placebo","Aspirin")
colnames(countTable) <- c("No heart attack","heart attack")
countTable

countTable["Placebo",]
countTable[,"No heart attack"]


bloodType <- c("A","B","AB","O","O","A","A","O","B","B")
summary(bloodType)
bloodType <- factor(bloodType)
bloodType
summary(bloodType)
levels(bloodType)
mode(bloodType)
class(bloodType)
gender <- c(1,1,2,2,1,2,2,2,1,2,2,1)
gender
summary(gender)
gender <- factor(gender)
gender
summary(gender)

gender <- c(1,1,2,2,1,2,2,2,1,2,2,1)
gender <- factor(gender,levels = c(1,2),labels = c("male","female")) # leveldmf label로 변환시킨다
gender
summary(gender)


###자료에 intake 라는 변수가있다.



intake.pre<- c(5260, 5470, 5640, 6180, 6390, 6515, 6805, 7515, 7515, 8230, 8770)
intake.post <-c(3910, 4220, 3885, 5160, 5645, 4680, 5265, 5975, 6790, 6990, 7335)

intake.race<- c('w', 'w', 'b', 'y', 'y', 'w', 'b', 'w', 'w', 'y', 'b')
intake <- data.frame(intake.pre,intake.post,intake.race)
intake
intake[,1]
intake$intake.pre

intake[1,]
colnames(intake) <- c("befere","after","race")
intake

intake <- data.frame(before = intake.pre,after = intake.post,race = intake.race)

head(intake)
tail(intake)


paste('w',1:11,sep = "")  # 이거 좋은듯...
paste('w',1:11,sep = "-")
rownames(intake)<- paste('w',1:11, sep ="")
intake['w4',]
summay(intake)

class(intake$race)
str(intake)


intake
