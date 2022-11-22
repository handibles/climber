
# your genomes go how fast??
wg <- data.frame( 
  g = c(32117,32245,32696,32834,33170,33435,34040,34510)-32116,
  t = c(18,21,34,37,44,50,63,74)-17
  )

plot(wg$g, wg$t, main = paste0("mean speed of ", round(mean(wg$g / wg$t),2) , " genomes per minute"), xlab = "n genomes", ylab = "n minutes")


(wgm <- lm(t ~ g, wg))
summary(wgm)
newd <- data.frame(g = 35000)
predict(wgm, newdata = newd)

