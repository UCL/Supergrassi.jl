# Supergrassi

[![Aqua QA](https://raw.githubusercontent.com/JuliaTesting/Aqua.jl/master/badge.svg)](https://github.com/JuliaTesting/Aqua.jl)
[![Build Status](https://github.com/UCL/Supergrassi.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/UCL/Supergrassi.jl/actions/workflows/CI.yml?query=branch%3Amain)

Multi-sector dynamic macroeconomics model with debt and default

## Read data into Julia

```julia
using Supergrassi

settings_path = create_filepath("config/settings.yml")
settings = read_settings(settings_path)
filepaths = check_file_availability(settings)
data = read_data(filepaths, settings)

```
