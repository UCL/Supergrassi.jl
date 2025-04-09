# settings_path = create_filepath("config/settings.yml")
# settings = read_settings(settings_path)
# filepaths = check_file_availability(settings)
# data = read_data(filepaths, settings)
clean = cleanup(data,settings["constants"]["data_year"])


println("Yaay")
# map_105_to_64 = create_map_105_to_64(data)

# println("Map 105 to 64: ", map_105_to_64)