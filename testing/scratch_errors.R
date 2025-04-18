# understand errors

# tryCatch(expression to try, error = something to do when there is an error,
# finally = something to always do?)
e <- simpleError("error")
tryCatch(e, error = function(e) e, finally = print("finally"))

