function Img=StackImage(Image)

for I=1:length(Image)
    Img(:,:,I)=Image{I};
end