# scatter plot : plot()

plot(2,1)
plot(c(2,5),c(1,10))

iris # 아이리스 데이터

plot(x=iris$Sepal.Length, y=iris$Sepal.Width)

plot(x=iris$Petal.Length, y=iris$Petal.Width)

cor(x=iris$Petal.Length, y=iris$Petal.Width)

# line plot : plot(typ='l','o')

plot(c(2,5),c(1,10),type = 'o')

pressure

plot(x=pressure$temperature,y=pressure$pressure,type='l')

# histogram: hist()

hist(iris$Sepal.Length)
hist(iris$Sepal.Width)

hist(iris$Sepal.Length, breaks = 20)

x <- matrix(rnorm(1000),nrow=100)
class(x)
hist(x)

# boxplot : boxplot()

boxplot(iris$Sepal.Length)
boxplot(iris$Sepal.Length,iris$Sepal.Width)

boxplot(iris[,1:4])

boxplot(iris$Sepal.Length~iris$Species)

boxplot(Petal.Length~ Species, data=iris)

# barplot : barplot()
head(BOD)

barplot(BOD$demand,names.arg=BOD$Time)
BOD

head(mtcars)
dim(mtcars)

cyl_freq = table(mtcars$cyl)
barplot(cyl_freq)

cyl_gear_freq = table(mtcars$cyl,mtcars$gear)
barplot(cyl_gear_freq)
barplot(cyl_gear_freq,beside = T)

# graphic parameter
barplot(cyl_gear_freq,beside=T
        ,xlab='Gears',ylab='Cylinder'
        ,ylim = c(0,20),las=1)

plot(x=pressure$temperature,y=pressure$pressure
      ,type='o',cex=0.5,col='red')

# par()

par(mfrow = c(2,2))
plot(x=pressure$temperature,y=pressure$pressure
     ,type='o',cex=0.5,col='blue')

graphics.off()
dev.off()

par(mar=c(10,10,10,10))

# ggplot
