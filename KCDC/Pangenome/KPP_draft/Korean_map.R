library(sf)
library(rnaturalearth)
library(ggplot2)
library(ggspatial)

install.packages(c("sf","rnaturalearth","rnaturalearthdata","ggspatial"))

world <- ne_countries(scale = "medium", returnclass = "sf")

# 2. 동북아 bbox 설정 (lon_min, lon_max, lat_min, lat_max)
bbox_neasia <- st_bbox(c(xmin = 110, xmax = 150, ymin = 5, ymax = 50), crs = st_crs(4326))

# 3. bbox로 자르기
neasia <- st_crop(world, bbox_neasia)

# 4. 기본 플롯
ggplot(neasia) +
  geom_sf(fill = "white", color = "black", size = 0.3) +
  geom_sf(data = st_intersection(neasia, neasia), fill = NA) + 
  #coord_sf(xlim = c(100,150), ylim = c(30,60), expand = FALSE) +
  annotation_scale(location = "bl") +
  annotation_north_arrow(location = "bl", which_north = "true", pad_x = unit(0.2,"cm")) +
  #labs(title = "동북아시아 (100E–150E, 30N–60N)") +
  theme(
    panel.background = element_blank(),    # 배경색 제거
    plot.background  = element_blank(),    # 플롯 전체 배경 제거
    panel.grid.major = element_blank(),    # 격자선 제거
    panel.grid.minor = element_blank(),
    axis.text        = element_blank(),    # 축 텍스트 제거
    axis.title       = element_blank(),    # 축 제목 제거
    axis.ticks       = element_blank(),    # 축 눈금 제거
    panel.border     = element_blank()
  )


skorea <- neasia[neasia$admin == "South Korea", ]

ggplot(neasia) +
  geom_sf(fill = "white", color = "#ecede6", size = 0.3) +        # 전체 지역 흰색
  geom_sf(data = skorea, fill = "#dde0d1", color = "#ecede6", size = 0.4) +  # 한국만 남색
  annotation_scale(location = "bl") +
  annotation_north_arrow(location = "bl",
                         which_north = "true",
                         pad_x = unit(0.2, "cm")) +
  theme(
    panel.background = element_blank(),
    plot.background  = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text        = element_blank(),
    axis.title       = element_blank(),
    axis.ticks       = element_blank(),
    panel.border     = element_blank()
  )
