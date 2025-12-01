sudo apt update
sudo apt-get install build-essential gdb
sudo apt-get install cmake
g++ --version
gdb --version


数量：100000 ✅
坐标：全部落在 [0,1]×[0,1] ✅
几何一致性：宽高、面积都自洽 ✅
覆盖率：≈ 0.736，和设定 0.7 很接近 ✅ 
大小分布：长尾、大小有明显差异 ✅
形状分布：多数接近正方形，也有长条矩形 ✅
重叠：平均每个矩形与约 4 个矩形相交 ✅


mkdir build
cd build
cmake ..
make -j
cd ..
./build/rect_sampler_experiments


git status
git add .
git commit -m "Initial commit"
git branch -M main
git push -u origin main


git fetch --all
git reset --hard origin/main
git pull

