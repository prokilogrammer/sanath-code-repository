
                       x <- seq(-3,3,by=.1)
y <- x^2
png(file="parabola.png")
           plot(x, y, xlab="x", ylab="y",
main="parabola y = x^2", type="l",
              col="blue")
dev.off() #turns off write-to-file device

