from readData import read_nodal_data
from obspy.core import UTCDateTime
from datetime import datetime

t1 = datetime.now()
event_loc = [10.85149,	-85.28347]
event_time=UTCDateTime("2023-07-23T00:18:22")
stations=[
    "453018446", 
   # "453016631", 
    #"453016442", 
    #"453015898",
   # "453014959",
    #"453014930",
    ]

p_arrivals=[
    "2023-07-23T00:18:23.01",
   #"2023-07-23T00:18:22",
   #"2023-07-23T00:18:21",
   #"2023-07-23T00:18:22",
   #"2023-07-23T00:18:22",
   # "2023-07-23T00:18:20",
    #"2023-07-23T00:18:19"
    ]

s_arrivals=[
    "2023-07-23T00:18:27.00",
   # "2023-07-23T00:18:25",
   # "2023-07-23T00:18:23",
   #"2023-07-23T00:18:26",
   # "2023-07-23T00:18:25",
   # "2023-07-23T00:18:25",
]

station_file = "nodos_coordenadas_select.txt"
for i, station in enumerate(stations): 
    read_nodal_data("data", station, station_file, phase_type='P&S', 
                    p_arrival=p_arrivals[i], s_arrival=s_arrivals[i], 
                    event_loc=event_loc, event_time=event_time)

t2 = datetime.now()
print(f"Execution time: {(t2-t1).total_seconds()} s" )