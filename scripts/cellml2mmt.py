import myokit.formats

model_dir = '../math_model/AP_model/'
i = myokit.formats.importer('cellml')
# m = i.model(model_dir + 'TT06.cellml')
# myokit.save_model(model_dir + 'TTP-2006.mmt', m)

# m = i.model(model_dir + 'Grd10.cellml')
# myokit.save_model(model_dir + 'Grd-2010.mmt', m)

m = i.model(model_dir + 'Tomek_model_epi.cellml')
myokit.save_model(model_dir + 'Tomek-2019.mmt', m)
