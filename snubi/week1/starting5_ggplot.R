# 변수 초기화
## rm(list=ls()) 

# 그래프 초기화
## graphics.off()
install.packages("ggplot2", dependencies=TRUE)


############ 3. Plot 함수 실습 

##### 들어가기 전에..
# y 값 할당
y <- 10:50 

# 할당된 값 확인
y        

# 그냥 그려보자
plot(y)


## lesson : 
# 1) 입력 data 값이 y 하나인데 2차원에 그려짐 - 선언되지 않은 x값은 y의 index값으로 잡혔다.
# 2) default가 선이 아니고 점이다.


# 실습에 배울 요소들 총합
y <- 1:5 
plot(y, main="This is\nmain", sub="This is sub", xlab="This is xlab", ylab="This is ylab", 
     type="o", lwd=2, col="blue", pch=24, bg="yellow", cex=1.5)


##### 1) plot type : type
par( mfrow=c(3, 3) )

#plot(y, type = 'p', main= "type P")

plot(y, type='p', main="type='p'" )
plot(y, type='l', main="type='l'" )
plot(y, type='b', main="type='b'" )
plot(y, type='c', main="type='c'" )
plot(y, type='o', main="type='o'" )
plot(y, type='h', main="type='h'" )
plot(y, type='s', main="type='s'" )
plot(y, type='S', main="type='S'" )
plot(y, type='n', main="type='n'" )

##### 2) line type : lty
graphics.off()
par( mfrow=c(2, 3) ) 

plot(y, type='o', lty=1, main=paste("lty=",1, sep='') )
plot(y, type='o', lty=2, main=paste("lty=",2, sep='') )
plot(y, type='o', lty=3, main=paste("lty=",3, sep='') )
plot(y, type='o', lty=4, main=paste("lty=",4, sep='') )
plot(y, type='o', lty=5, main=paste("lty=",5, sep='') )
plot(y, type='o', lty=6, main=paste("lty=",6, sep='') )


for(i in 1:6)
  plot(y, type='o', lty=i, main=paste("lty=",i, sep='') )

##### par(mfrow =), par(mfcol =)
graphics.off()
par( mfrow=c(2, 3) ) 

plot(y, type='o', lty=1, main=paste("lty=",1, sep='') )
plot(y, type='o', lty=2, main=paste("lty=",2, sep='') )
plot(y, type='o', lty=3, main=paste("lty=",3, sep='') )
plot(y, type='o', lty=4, main=paste("lty=",4, sep='') )
plot(y, type='o', lty=5, main=paste("lty=",5, sep='') )
plot(y, type='o', lty=6, main=paste("lty=",6, sep='') )


graphics.off()
par( mfcol=c(2, 3) ) 

plot(y, type='o', lty=1, main=paste("lty=",1, sep='') )
plot(y, type='o', lty=2, main=paste("lty=",2, sep='') )
plot(y, type='o', lty=3, main=paste("lty=",3, sep='') )
plot(y, type='o', lty=4, main=paste("lty=",4, sep='') )
plot(y, type='o', lty=5, main=paste("lty=",5, sep='') )
plot(y, type='o', lty=6, main=paste("lty=",6, sep='') )


##### 4) plot character/symbol : pch 
par( mfrow=c(3, 3) )
for ( i in 1:9 ) plot(y, type='o', pch=i, main=paste('pch=', i, sep='') )

par(mfrow = c(1,2)) 
plot(y, type="o", col="blue")
plot(y, type="o", pch=23, col="red3", bg="slateblue3")

plot(y, type="o", pch=23, col="red", bg="blue3", col.axis = "cyan", col.main="blue2", main = "strange")

##### 번외 : color 지정 방법 4가지 
par(mfrow = c(1,4)) 
# 1) 숫자 (0-8까지만 가능)
plot(y, col=6, cex = 2, pch = 19)
# 2) color name
plot(y, col="red", cex = 2, pch = 19)
# 3) 16진법 표기
plot(y, col="#a8259f", cex = 2, pch = 19)
# 4) RGB 표기
rgb <- rgb(168,37,159, maxColorValue=255)
plot(y, col=rgb, cex = 2, pch = 19)

# col에 변수를 할당해서 알록달록 그리기 (반복 색할당)
rgb <- c("red","green","blue")
plot(y, col=rgb, cex = 2, pch = 19)


###############################################################
##################     Low Level Plots         ################
###############################################################

##### 5. Lines, axis, box, title, legend, matplot 함수

# par( mfrow=c(1, 1) )

# Dataset 만들기
cars <- data.frame( 
  standard=c(1, 3, 6, 4, 9),
  truck=c(2, 5, 4, 5, 12),
  suv=c(4, 4, 6, 6, 16) )

rownames(cars) <- c("Mon", "Tue", "Wed", "Thu", "Fri")
cars

# plot(cars) 
# 차종별로 각각 그려보기
par( mfrow=c(1, 3) )

plot( cars$standard, type="o")
plot( cars$truck, type="o")
plot( cars$suv, type="o")

# 빈 plot부터 요소를 하나씩 그려보자
par( mfrow=c(1, 1) )

plot( cars$standard, type="o", col="blue", ylim=c(0, max_y), axes=FALSE, ann=FALSE )

plot( cars$standard, type="o", col="blue", ylim=c(0, max_y))
plot( cars$standard, type="o", col="blue", ylim=c(0, max_y), axes=FALSE)
plot( cars$standard, type="o", col="blue", ylim=c(0, max_y), axes=FALSE, ann=FALSE )

# 축을 그리자 : axis(), box() 함수
axis(1, at=1:5, lab=rownames(cars) )    # axis 1 x축 / 2 y축



max_y <- max(cars)
max_y

axis(2, at=seq(0, max_y, by=4), las=1 ) # las = 1 가로 / las = 2 세로
box()


# 선을 그리자 : lines() 함수
lines(cars$truck, type="o", pch=22, lty=2, col="red")
lines(cars$suv, type="o", pch=23, lty=3, col="green")

#범례 그리기 : legend() 함수
title(main="Car Rental", xlab="Weekday", ylab="The Number of Cars")
legend( "topleft", colnames(cars), col=c("blue", "red", "green"), pch=21:23, lty=1:3 )

##### 기타 범례 위치 지정
plot(1:10, type ="n", main = "LEGEND")

legend( "bottom", colnames(cars), col=c("blue", "red", "green"), pch=21:23, lty=1:3, title="bottom" )
legend( "bottomleft", colnames(cars), col=c("blue", "red", "green"), pch=21:23, lty=1:3, title="bottomleft" )
legend( "bottomright", colnames(cars), col=c("blue", "red", "green"), pch=21:23, lty=1:3, title="bottomright" )
legend( "left", colnames(cars), col=c("blue", "red", "green"), pch=21:23, lty=1:3, title = "left" )
legend( "right", colnames(cars), col=c("blue", "red", "green"), pch=21:23, lty=1:3, title = "right" )
legend( "top", colnames(cars), col=c("blue", "red", "green"), pch=21:23, lty=1:3, title = "top" )
legend( "topleft", colnames(cars), col=c("blue", "red", "green"), pch=21:23, lty=1:3, title = "topleft" )
legend( "topright", colnames(cars), col=c("blue", "red", "green"), pch=21:23, lty=1:3, title = "topright" )
legend( "center", colnames(cars), col=c("blue", "red", "green"), pch=21:23, lty=1:3, title = "center" )

### 좌표 입력 
legend(7,8, legend=colnames(cars), col=c("blue", "red", "green"), pch=21:23, lty=1:3, title="(x,y)" )

# 마우스로 위치 찍어서 그리기 
legend(locator(1), legend="Locator", fill = 1 ) # 상자 크기 text에 맞추기
legend(locator(2), legend="Locator", fill = 1 ) # 상자 크기 지정


#########################################################################
############# attach(), detach(), hist(), density(), boxplot() ##########
#########################################################################

## Dataset birthwt 사용 
library(MASS)

head(birthwt)
dim(birthwt)
summary(birthwt)

?birthwt

factor(birthwt$race)

birthwt$race <- factor( birthwt$race, levels=c(1, 2, 3), labels=c("white", "black", "other") )
levels(birthwt$race)

## attach(), detach()
head(birthwt$bwt)
head(bwt)

attach(birthwt)
head(bwt)

detach(birthwt)

## hist() : age변수를 탐색해보자     
attach(birthwt)
summary(birthwt)
dim(birthwt)
?birthwt

head(birthwt)
plot(birthwt) 

plot(age, main = "Age")

hist(age)

# 교재의 예제를 통해 파라메터를 익혀보자
par( mfrow=c(2, 3) )
hist(age)
hist(age, freq=FALSE)

hist(age, breaks=c(10, 15, 17, 19, 21, 24, 25, 30, 35, 40, 45) )
hist(age, labels=c("1~14", "15~19", "20~24", "25~29", "30~34", "35~39", "40~44") )
hist(age, labels=c("1~14", "15~19", "20~24", "25~29", "30~34", "35~39", "40~44"), col="red")
hist(age, labels=c("1~14", "15~19", "20~24", "25~29", "30~34", "35~39", "40~44"), col="red", 
     density=10)

# 데이터가 잘 파악되는 그림을 찾자
par(mfrow=c(1,1))
hist(age, breaks=c(10,20,30,40,50), labels=c("10대","20대","30대", "40대"))
hist(age, breaks=c(10,15,20,25,30,35,40,45,50), col = 2:9, density = 40)

## density 함수: bwt의 분포를 확인해보자
par( mfrow=c(1,1))
plot( density(bwt))

plot( density(bwt[race=='black']), col="black" )
lines( density(bwt[race=='white']), col="red" )
lines( density(bwt[race=='other']), col="blue" )
legend( "topright", legend=levels(race), fill=c("red", "black", "blue") )

## boxplot : 인종별 bwt의 차이? 
head( birthwt[, c("bwt", "race")] )

par( mfrow=c(1, 2) )
boxplot(bwt)
boxplot(bwt ~ race)

boxplot(bwt~smoke, col=rainbow(2), xlab="Smoke", ylab="Birth Weight")
boxplot(bwt~race, col=rainbow(3), xlab="Race", ylab="Birth Weight")
detach(birthwt)

## boxplot
par( mfrow=c(1, 1) )
barplot( c(1, 2, 3, 4), names.arg=c("A", "B", "C", "D") )

head(iris)
str(iris)

attach(iris)

iris.mean <- aggregate(iris[,1:4], iris["Species"], mean)
iris.mean
iris.sd <- aggregate(iris[,1:4], iris["Species"], sd)
iris.sd
iris.sd.upper <- iris.mean[,2:5] + iris.sd[,2:5]
iris.sd.lower <- iris.mean[,2:5] - iris.sd[,2:5]

par( mfrow=c(1, 1) )
b <- barplot( as.matrix(iris.mean[,2:5]), beside=TRUE, col=c("red","blue","purple"), ylim=c(0,8) )
arrows(b, as.matrix(iris.sd.upper), b, as.matrix(iris.sd.lower), angle=90, length=0.05, code=3)
legend( "topright", legend=iris.mean$Species, fill=c("red","blue","purple") )

###############################################################
########    8. ggplot / introduction & install      ###########
###############################################################


install.packages("ggplot2", dependencies=TRUE)
library(ggplot2)

graph <- ggplot( iris, aes(Sepal.Length, Sepal.Width) )
graph + geom_point()

graph + geom_point( aes(color=Species, shape=Species), size=3 )
graph + geom_point( aes(color=Species, size=Sepal.Length/Sepal.Width) )

graph <- ggplot( iris, aes(Species, Sepal.Length, fill=Species) )
graph + geom_boxplot() 
graph + geom_boxplot() + geom_jitter()

