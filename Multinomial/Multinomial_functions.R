#############################################################
# MultinomialLasso
#############################################################

trainMnLasso <- function(Expression, Y) {
  Salida <- cv.glmnet(Expression,Y, family = "multinomial", 
                      alpha = 1, type.multinomial = "grouped")
  return(Salida)
}

predictMnLasso <- function(Salida, Expression) {
  Prediccion <- drop(predict(Salida,Expression, s="lambda.min"))
  TratamientoMnLasso <- apply(Prediccion, 1, which.max)
  TratamientoMnLasso <- factor(TratamientoMnLasso, levels = 1:ncol(C))
  return(TratamientoMnLasso)
}
