# Temperature through time

Using the Kawamura ice core data from [here](https://www.ncei.noaa.gov/pub/data/paleo/icecore/antarctica/domefuji/df-tsite-340ka-dfo2006.txt) 

Additional data on the dataset [here](https://www.ncei.noaa.gov/access/metadata/landing-page/bin/iso?id=noaa-icecore-6076). 

## Assign Two Hot & Cold Periods

Quantitatively assign t0 and t1 hot and cold periods. 

For example: 

```
# Hot times, recent
hot_0 <- icecore %>%  filter(time < 5e4) %>% slice_max(deltaT,n=20) %>%
  summarize(time_low = 0, time_high = max(time), temp = mean(deltaT), duration = time_high - time_low)

# Cold times, recent
cold_0 <- icecore %>%  filter(time < 5e4) %>% slice_min(deltaT,n=20) %>%
  summarize(time_low = min(time), time_high = max(time), temp = mean(deltaT), duration = time_high - time_low)

# Hot times, mid
hot_1 <- icecore %>%  filter(time < 1.75e5 & time > 5e4) %>% slice_max(deltaT,n=20) %>%
  summarize(time_low = min(time), time_high = max(time), temp = mean(deltaT), duration = time_high - time_low)

# Cold times, mid
cold_1 <- icecore %>%  filter(time < 1e5 & time > 5e4) %>% slice_min(deltaT,n=20) %>%
  summarize(time_low = min(time), time_high = max(time), temp = mean(deltaT), duration = time_high - time_low)

```

But, we want our time periods to be a similar duration, we will select 10Ka for simplicity. I find the midpoint of the periods, and then add +/- 5000 so that all periods are the same duration (10K). Also ensure that the contemporary period starts at 0 and goes to 10Ka. 

| Start  | End    | Duration | Stage | Period | Mean Temperature |
| ------ | ------ | -------- | ----- | ------ | ---------------- |
| 0      | 10000  | 10000    | hot   | t0     | 0.005            |
| 23750  | 33750  | 10000    | cold  | t0     | -7.6             |
| 121500 | 131500 | 10000    | hot   | t1     | 1.94             |
| 56750  | 66750  | 10000    | cold  | t1     | -6.46            |


Which gives us:

![Periods](/figures/WarmCold_Periods.png)

