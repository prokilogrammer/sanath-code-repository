A = read.table("input")
jpeg("isotope-profile.jpg")
plot(A[1],A[2] ,col="blue1", main="Isotope Profile for SLAMMER", xlab='Isotope Peak #', ylab='Probability', labels=TRUE)
dev.off()

