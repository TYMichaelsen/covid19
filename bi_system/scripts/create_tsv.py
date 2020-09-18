import csv

def save_csv(path, data, fields):
      with open(path, "w") as f:
            writer = csv.DictWriter(f, fieldnames=fields, delimiter='\t')
            writer.writeheader()
            for data_obj in data:
                  writer.writerow(data_obj)


fields = ['code', 'name', 'latitude', 'longitude', 'nuts2_region']
path = './bi_system/stable_dims/nuts3_regions.tsv'
data = [
    {'code':'DK011', 'name':'Byen København', 'latitude':'55.61111592', 'longitude':'12.59840382', 'nuts2_region':'DK01'},
    {'code':'DK012', 'name':'Københavns omegn', 'latitude':'55.66098406', 'longitude':'12.3885576', 'nuts2_region':'DK01'},
    {'code':'DK013', 'name':'Nordsjælland', 'latitude':'55.97239387', 'longitude':'12.27937981', 'nuts2_region':'DK01'},
    {'code':'DK014', 'name':'Bornholm', 'latitude':'55.12789092', 'longitude':'14.8836754', 'nuts2_region':'DK01'},
    {'code':'DK021', 'name':'Østsjælland', 'latitude':'55.54415661', 'longitude':'12.06527713', 'nuts2_region':'DK02'},
    {'code':'DK022', 'name':'Vest- og Sydsjælland', 'latitude':'55.41305103', 'longitude':'11.5746375', 'nuts2_region':'DK02'},
    {'code':'DK031', 'name':'Fyn', 'latitude':'55.23531386', 'longitude':'10.45714094', 'nuts2_region':'DK03'},
    {'code':'DK032', 'name':'Sydjylland', 'latitude':'55.56316692', 'longitude':'9.02846054', 'nuts2_region':'DK03'},
    {'code':'DK041', 'name':'Vestjylland', 'latitude':'56.18780543', 'longitude':'8.73661729', 'nuts2_region':'DK04'},
    {'code':'DK042', 'name':'Østjylland', 'latitude':'56.1071088', 'longitude':'9.81621671', 'nuts2_region':'DK04'},
    {'code':'DK051', 'name':'Nordjylland', 'latitude':'57.30715928', 'longitude':'10.11282907', 'nuts2_region':'DK05'}
]
save_csv(path, data, fields)

fields2 = ['code', 'name', 'latitude', 'longitude']
path2 = './bi_system/stable_dims/nuts2_regions.tsv'
data2 = [
    {'code':'DK01', 'name':'Region_Hovedstaden', 'latitude':'55.97238879', 'longitude':'12.27937393'},
    {'code':'DK02', 'name':'Region_Sjælland', 'latitude':'55.43978846', 'longitude':'11.6213381'},
    {'code':'DK03', 'name':'Region_Syddanmark', 'latitude':'55.56316692', 'longitude':'9.02846054'},
    {'code':'DK04', 'name':'Region_Midtjylland', 'latitude':'56.23399042', 'longitude':'9.60503157'},
    {'code':'DK05', 'name':'Region_Nordjylland', 'latitude':'57.30715928', 'longitude':'10.11282907'}
]
save_csv(path2, data2, fields2)