# Подключим необходимые библиотеки
library("titanic")                                               # пакет, содержащий датафрейм с
                                                                 # информацией о пассажирах Титаника
library("glmx")                                                  # пакет, позволяющий оценивать пробит
                                                                 # модель с гетероскдестичной 
                                                                 # случайной ошибкой
library("lmtest")                                                # дополнительные тесты

library("crch")
library("hpa")

library("numDeriv")                                              # численное дифференцирование
library("sampleSelection")                                       # встроенный датафрейм для
                                                                 # выполнения домашнего задания
library("GJRM")                                                  # система бинарных уравнений

# Отключим scientific notation
options(scipen = 999)

# Подключим встроенный датафрейм, содержащий информацию
# о характеристиках американских женщин и их трудовом
# статусе в 1975-1976 годах
data(Mroz87)
# Для тех, у кого не считывается встроенный датафрейм
library("readxl")                                                                          # предварительно установите пакет
Mroz87 <- read_excel("G:\\Преподавание\\Микроэконометрика\\Домашнее задание\\Mroz87.xlsx") # укажите путь к файлу                                                                                          # на своем компьютере
h <- na.omit(as.data.frame(Mroz87))
# Краткое описание:
help(Mroz87)
# lfp - бинарная переменная, принимающая значение 1
#       если женщина работает и 0 - в противном случае
# hours - количество проработанных часов за год
# kids5 - количество детей младше 6 лет
# kids618 - количество несовершеннолетних детей
#           старше пяти лет
# age - возраст женщины
# educ - число лет, потраченных женщиной на
#        получение образования
# wage - почасовая зарплата женщины
# hushrs - количество часов, проработанных мужем женщины
# husage - возраст мужа женщины
# huseduc - число лет, потраченных муженм женщины на
#           получение образования
# huswage - зарплата мужа женщины
# faminc - доход семьи женщины
# mtr - налоговая нагрузка на женщину
# fatheduc - число лет, потраченных отцом женщины на
#            получение образования
# motheduc - число лет, потраченных отцом женщины на
#            получение образования
# unem - безработица в регионе проживания женщины
# city - бинарная переменная, принимающая значение 1
#        если женщина живет в городе и 0 - иначе
# exper - рабочий стаж женщины в годах
# nwifeinc - доход семьи женщины за вычетом ее дохода
# wifecoll - бинарная переменная, принимающая значение 1
#            если женщина посещала колледж и 0 - иначе
# huscoll - бинарная переменная, принимающая значение 1
#           если муж женщины посещал колледж и 0 - иначе


# 2. МОДЕЛЬ ТОБИНА

# 2.1. - 2.3.

tr_left <- 0                                                     # наблюдаем з/п только для тех,
                                                                 # у кого она есть (= положительная)
tr_right <- Inf                                                  # справа цензурирование отсутствует

model_tobit <- crch(wage ~ age + I(age^2) + city + unem,   
                    data = h,                        
                    left = tr_left,                         
                    right = tr_right)    

summary(model_tobit)


# 2.4.
# Индивид с произвольными характеристиками

individ <- data.frame("age" = 35,
                     "city" = 1,
                     "unem" = 9)

# a) E(y*)

wage_star_est <- predict(model_tobit, 
                         newdata = individ) 

# б) Е(у)

est_tobit <- coef(model_tobit)                                   # оценки коэффициентов модели
coef_tobit <- est_tobit[-length(est_tobit)]    
sigma_tobit <- exp(est_tobit[length(est_tobit)])                 # оценка стандартного отклонения

epsilon_E <- truncatedNormalMoment(k = 1,                        # условное математическое ожидание      
                                   x_lower = tr_left -           # случайной ошибки
                                             wage_star_est, 
                                   x_upper = tr_right - 
                                             wage_star_est,
                                   mean = 0, sd = sigma_tobit)

prob_est_left <- predict(model_tobit,                            # оценка вероятности усечения снизу
                         at = tr_left,                  
                         type = "probability", 
                         newdata = individ)    

wage_est <- (wage_star_est + epsilon_E) * (1 - prob_est_left) +  # E(y|y* > lower) * P(y* > lower)
            tr_left *  prob_est_left                             # lower * P(y* < lower)           
                 
# в) P(y* > lower

prob_est_right <- 1 - prob_est_left


## 2.5.
# а) предельный эффект на E(y*)

ME_wage_star_unem <- coef_tobit["unem"]

# б) предельный эффект на E(y)

ME_wage_unem <- coef_tobit["unem"] * prob_est_right

# в) предельный эффект на P(y* > lower)

ME_prob_right_unem <- coef_tobit["unem"] *           
                        dnorm(wage_star_est - tr_left,      
                              sd = sigma_tobit)  


## 2.6.
# а) предельный эффект на E(y*)

ME_wage_star_age <- coef_tobit["age"] + 
                    2 * individ$age * coef_tobit["I(age^2)"]

# б) предельный эффект на E(y)

ME_wage_age <- (coef_tobit["age"] + 
                2 * individ$age * coef_tobit["I(age^2)"]) *
                prob_est_right
# в) предельный эффект на P(y* > lower)

ME_prob_right_age <- (coef_tobit["age"] + 
                     2 * individ$age * coef_tobit["I(age^2)"]) *           
                     dnorm(wage_star_est - tr_left,      
                           sd = sigma_tobit)  

# 3. МОДЕЛЬ ХЕКМАНА 
## 3.1.

model_heckman_mle <- selection(
                       selection = lfp ~ huswage + mtr + unem,   # уравнение отбора
                       outcome = wage ~ age + I(age^2) +         # основное уравнение
                                        city + unem,                                
                       data = h,                                       
                       method = "ml")  

summary(model_heckman_mle)      

coef_mle <- coef(model_heckman_mle, part = "outcome")             # оценки коэффициентов
rho_mle <- model_heckman_mle$estimate["rho"]                      # оценка корреляции между
                                                                  # случайными ошибками


## 3.3.

model_heckman_2st <- selection(                              
                        selection = lfp ~ huswage + mtr + unem,                     
                        outcome = wage ~ age + I(age^2) +         
                                         city + unem,                  
                        data = h,                                          
                        method = "2step")    

summary(model_heckman_2st)                                 

coef_2st <- coef(model_heckman_2st, part = "outcome")             # оценки коэффициентов
coef_2st <- coef_2st[-length(coef_2st)]              
rho_2st <- model_heckman_2st$rho                                  # оценка корреляции между
                                                                  # случайными ошибками

# Сравнение методов
# регрессионные коэффициенты

coef_compare <- data.frame("Heckman MLE" = coef_mle,                 
                           "Heckman 2step" = coef_2st)   

# корреляция случайных ошибок

corr_compare <- data.frame("Heckman MLE" = rho_mle,                  
                           "Heckman 2step" = rho_2st)                


## 3.5.

individ$lfp <- 1                                                   # добавление недостающих 
individ$wage <- 10                                                 # характеристик
individ$huswage <- 13
individ$mtr <- 0.55

# a) E(y*|z=1) и E(y*|z=0)

wage_cond <- predict(model_heckman_mle, 
                     newdata = individ, 
                     part = "outcome",                   
                     type = "conditional")               
wage_cond[1]                                                      # E(y*|z = 0)
wage_cond[2]                                                      # E(y*|z = 1)

# б) предельный эффект на E(y*|z=1) и E(y*|z=0)

eps <- 1e-8
individ_eps <- individ
individ_eps$unem <- individ$unem + eps
wage_cond_eps <- predict(model_heckman_mle, 
                         newdata = individ_eps, 
                         part = "outcome",
                         type = "conditional")
unem_ME_0 <- (wage_cond_eps[1] - wage_cond[1]) / eps            # E(y*|z = 0)
unem_ME_1 <- (wage_cond_eps[2] - wage_cond[2]) / eps            # E(y*|z = 1)









