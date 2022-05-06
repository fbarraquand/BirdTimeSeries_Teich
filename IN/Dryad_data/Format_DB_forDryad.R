### FB 06/05/2022 - removing remarks and first column
DB<-read.csv(file="../DBWithMonthlyPhotoTeich_completed.csv",header=TRUE,sep=",",dec=".")
DB$X<-NULL
DB$Commentaire<-NULL
DB$Remarque_privee<-NULL
label = "DBWithMonthlyPhotoTeich" # new name of the database
write.csv(DB,file = paste(label,".csv",sep=""))
