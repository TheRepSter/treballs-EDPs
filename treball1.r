# Pràctica 1
# Eduardo Pérez Motato 1709992
# Fèlix Sáiz von Fraunberg 1620854

####### Exercici 1 #######
# a) Trobar la solució u(a,t) en termes de b(t) i u0(a)
# 
# L'equació és: ∂_t u + ∂_a u = -c(a) u, amb condicions u(0,t)=b(t), u(a,0)=u0(a)
# Utilitzant el mètode de les característiques:
# - Si t < a: la característica ve de la condició inicial, u(a,t) = u0(a-t) * exp(-∫_{a-t}^a c(s) ds)
# - Si t >= a: la característica ve de la condició de contorn, u(a,t) = b(t-a) * exp(-∫_0^a c(s) ds)
#
# b) Per t > 1, com que a <= 1 < t, sempre es dóna t >= a, per tant u(a,t) = b(t-a) * exp(-∫_0^a c(s) ds), que només depèn de b(t).
#
# c) Si b(t) = b constant, per t > 1, u(a,t) = b * exp(-∫_0^a c(s) ds) = ū(a). 
# Comprovem que ū(a) satisfà el problema estacionari:
#   u'(a) = -c(a) u(a), u(0)=b.
# Derivant: ū'(a) = b * exp(-∫_0^a c(s) ds) * (-c(a)) = -c(a) ū(a)
# i ū(0) = b * exp(0) = b.
# Per tant, és solució.

####### Exercici 2 #######

# Dades d'entrada
u0 <- function(a) { a }
cf <- function(a) { 1 / (1 - a) }  # Mortalitat c(a)
b <- function(t) { 1 }              # Naixements constants
M <- 500
da <- 1 / M
mu <- 1/2
dt <- mu * da
t_final <- 1

# Malla espacial
av <- (0:(M-1)) * da   # Punts d'edat: a_0=0, a_1=da, ..., a_{M-1}=1-da

# Inicialització del vector U (condició inicial)
U <- u0(av)

# Nombre de passos de temps
n_steps <- round(t_final / dt)

# Bucle per als passos de temps (esquema explícit)
for (n in 1:n_steps) {
  U_new <- numeric(M)
  # Condició de contorn per a m=0 (índex 1 en R)
  U_new[1] <- b(n * dt)
  # Per a m=1,...,M-1 (índexs 2 a M en R)
  for (i in 2:M) {
    c_val <- cf(av[i])
    U_new[i] <- (1 - dt/da - dt * c_val) * U[i] + (dt/da) * U[i-1]
  }
  U <- U_new
}

# Solució numèrica a t=1
U_numeric <- U

# Solució analítica a t=1: u(a,1) = 1 - a (com demostrat a l'exercici 1)
U_analytic <- 1 - av

####### Exercici 3 #######
plot(av, U_analytic, type = "l", col = "black", lwd = 2, 
     xlab = "Edat (a)", ylab = "Densitat u(a,1)", 
     main = "Comparació de la solució numèrica i analítica a t=1")
lines(av, U_numeric, col = "red", lty = 2, lwd = 2)
legend("topright", legend = c("Analítica", "Numèrica"), 
       col = c("black", "red"), lty = c(1,2), lwd = 2)

####### Exercici 4 #######

# Dades amb Δt = 2Δa
M <- 500
da <- 1 / M
dt <- 2 * da   # μ = 2
t_final <- 1
n_steps <- round(t_final / dt)

# Malla espacial
av <- (0:(M-1)) * da

# Inicialització
U <- u0(av)

# Trobar l'índex per a a=0.5
index_a05 <- which.min(abs(av - 0.5))  # Ha de ser 251

# Calcular k per mostrejar cada 0.1 unitats de temps
k <- floor(t_final / (10 * dt))   # k = 25

# Vector per emmagatzemar els valors a a=0.5 en temps seleccionats
temps <- seq(0, 1, by = 0.1)
valors_U <- numeric(length(temps))

# Guardar el valor inicial
valors_U[1] <- U[index_a05]

# Bucle de passos de temps
for (n in 1:n_steps) {
  U_new <- numeric(M)
  U_new[1] <- b(n * dt)
  for (i in 2:M) {
    c_val <- cf(av[i])
    U_new[i] <- (1 - dt/da - dt * c_val) * U[i] + (dt/da) * U[i-1]
  }
  U <- U_new
  # Guardar si n és múltiple de k
  if (n %% k == 0) {
    idx <- n / k + 1   # Índex per a valors_U (començant a t=0.1, 0.2, ...)
    if (idx <= length(valors_U)) {
      valors_U[idx] <- U[index_a05]
    }
  }
}

taula <- data.frame(temps = temps, u_05 = valors_U)
print(taula)

####### Exercici 5 #######

# Funció per resoldre el model per a un r donat
model_r <- function(r, M, t_final) {
  da <- 1 / M
  # Triem Δt = 0.5Δa per estabilitat
  dt <- 0.5 * da
  n_steps <- round(t_final / dt)
  av <- (0:(M-1)) * da
  U <- u0(av)   # Condició inicial u0(a)=a
  
  for (n in 1:n_steps) {
    # Calcul de la població total P_n
    P_n <- sum(U) * da
    # Condició de contorn dependent de P_n
    b_n <- r * P_n / (1 + P_n)
    U_new <- numeric(M)
    U_new[1] <- b_n
    for (i in 2:M) {
      c_val <- cf(av[i])
      U_new[i] <- (1 - dt/da - dt * c_val) * U[i] + (dt/da) * U[i-1]
    }
    U <- U_new
  }
  # Retorna la població total a t_final
  return(sum(U) * da)
}

# Paràmetres per a l'exercici 5
M <- 100
t_final <- 50
vals_r <- seq(0, 4, by = 0.1)
P_final <- sapply(vals_r, model_r, M = M, t_final = t_final)

# Gràfic de P(50) en funció de r
plot(vals_r, P_final, type = "o", pch = 16, cex = 0.5, 
     xlab = "Paràmetre r", ylab = "Població total P(50)",
     main = "Població total en funció de la fertilitat r")

# Comentari: Es pot veure que per a r petit, P(50) és zero, i per a r gran, P(50) és positiva.
# La condició teòrica per a la persistència és ∫_0^1 r e^{-∫_0^a c(α)dα} da > 1.
# En el nostre cas, c(a)=1/(1-a), per tant ∫_0^a c(α)dα = -ln(1-a), i e^{-∫_0^a c(α)dα} = 1-a.
# Aleshores ∫_0^1 r (1-a) da = r [a - a^2/2]_0^1 = r/2.
# La condició és r/2 > 1, és a dir, r > 2. 
# A la gràfica s'observa que per r > 2, P(50) és positiva, consistent amb la teoria.