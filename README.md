
#INSTALLATION

#BUILD
docker build -t liquid_growth_analysis .

#FILE TREE

#EXECUTE
docker run --name 2303_21 -v /Users/groot/Desktop/growthcurving/:/workdir -it liquid_growth_analysis 

#END
docker rm $(docker ps -qa)