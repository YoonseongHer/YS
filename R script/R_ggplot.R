ggplot(data=iris,aes(x=Petal.Length,y=Petal.Width)
       )+geom_point(size=3,color='red')

# geom_line
ggplot(data=pressure,
       aes(x=temperature,y=pressure)
       ) +
  geom_line()+geom_point()
