dir = "./datas"
file = readdir(dir)
for f in file
    rm(joinpath(dir,f))
end

dir = "./datasCyl"
file = readdir(dir)
for f in file
    rm(joinpath(dir,f))
end