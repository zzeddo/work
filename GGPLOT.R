#------------------------------------------------
# 1. 소개
library(ggplot2)

#------------------------------------------------
# 2. 데이터셋
# ggplot2 diamonds dataset(54,000 diamonds values)
# columns are : carat, cut, color, clarity, depth, price, x, y
# depth = z depth / z * 100
# Table = table width / x * 100

diamonds[1:3,]
#  carat     cut color clarity depth table price    x    y    z
#1  0.23   Ideal     E     SI2  61.5    55   326 3.95 3.98 2.43
#2  0.21 Premium     E     SI1  59.8    61   326 3.89 3.84 2.31
#3  0.23    Good     E     VS1  56.9    65   327 4.05 4.07 2.31

# Random sample of 100 diamonds.
set.seed(1410) # Make the sample reproducible
dsmall <- diamonds[sample(nrow(diamonds), 100), ]

#------------------------------------------------
# 3. 기본 사용
ggplot(diamonds, aes(x=carat, y=price))+ geom_point()

# 이상치(Outliers)로 인한 왜곡을 줄이기 위하여 log(price)와 log(carat)을 통하여 재확인
ggplot(diamonds, aes(x=log(carat), y=log(price)))+ geom_point()

# 다이아몬드의 3차원 좌표(x, y, z)에 대한 면적(Valume : x*y*z)와 무게(carat)과의 관계 도식화
ggplot(diamonds, aes(x=carat, y=x*y*z))+ geom_point()

#------------------------------------------------
# 4. 색깔(color), 크기(size), 모양(shpe)와 다른 가시적(aesthetic) 속성


ggplot(diamonds, aes(x=carat, y=price, color=cut))+ geom_point()
