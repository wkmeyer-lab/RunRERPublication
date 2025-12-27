fahrenheitToKevlin = function(input) {
  #converts input of fahrenhiet to Kevlin
  stopifnot(is.numeric(input))
  celsius = ((input -32)* (5/9))
  kelvin = celsius + 273.15
  return(kelvin)
}

fahrenheitToCelsius = function(input) {
  #converts input of fahrenhiet to Celcius
  stopifnot(is.numeric(input))
  celsius = ((input -32)* (5/9))
  kelvin = celsius + 273.15
  return(celsius)
}

kelvinToCelsius = function(input){
  #converts input of Kevlin to Celsius
  stopifnot(is.numeric(input))
  celsius = input - 273.15
  return(celsius)
}

celsiusToKevlin = function(input){
  #converts input of Celsius to Kelvin
  stopifnot(is.numeric(input))
  kelvin = input + 273.15
  return(kelvin)
}
celsiusToFahernhiet = function(input){
  #converts input of Celsius to Fahrenhiet
  stopifnot(is.numeric(input))
  fahernheit= input * (9/5) +32
}
kelvinToFahernhiet = function(input){
  #converts input of Kelvin to Fahrenhiet
  stopifnot(is.numeric(input))
  celsius = input - 273.15
  fahernheit = celsius * (9/5) +32
}