function varargout = Image_Fusion(varargin)
% IMAGE_FUSION MATLAB code for Image_Fusion.fig
%      IMAGE_FUSION, by itself, creates a new IMAGE_FUSION or raises the existing
%      singleton*.
%
%      H = IMAGE_FUSION returns the handle to a new IMAGE_FUSION or the handle to
%      the existing singleton*.
%
%      IMAGE_FUSION('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in IMAGE_FUSION.M with the given input arguments.
%
%      IMAGE_FUSION('Property','Value',...) creates a new IMAGE_FUSION or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before Image_Fusion_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to Image_Fusion_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help Image_Fusion

% Last Modified by GUIDE v2.5 13-May-2020 11:15:44

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @Image_Fusion_OpeningFcn, ...
                   'gui_OutputFcn',  @Image_Fusion_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before Image_Fusion is made visible.
function Image_Fusion_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to Image_Fusion (see VARARGIN)

% Choose default command line output for Image_Fusion
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes Image_Fusion wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = Image_Fusion_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in pushbutton1.
function pushbutton1_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[filename1,PathName1] = uigetfile({'*.BMP';'*.bmp';'*.tif';'*.jpg';'*.png'},'C:\Users\韩毅\Documents\MATLAB\multi-focus');
X1 = [PathName1 filename1];
if PathName1 ~=0
    OriginImage1 = imread(X1);
    handles.OrginImage1=OriginImage1;
    guidata(hObject,handles);
    axes(handles.axes1);
    imshow(OriginImage1);
end

% --- Executes on button press in pushbutton2.
function pushbutton2_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[filename2,PathName2] = uigetfile({'*.BMP';'*.bmp';'*.tif';'*.jpg';'*.png'},'C:\Users\韩毅\Documents\MATLAB\multi-focus');
X2 = [PathName2 filename2];
if PathName2 ~=0
    OriginImage2 = imread(X2);
    handles.OrginImage2=OriginImage2;
    guidata(hObject,handles);
    axes(handles.axes2);
    imshow(OriginImage2);
end

% --- Executes on button press in pushbutton3.
function pushbutton3_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%图像融合
OriginImage1=handles.OrginImage1;
OriginImage2=handles.OrginImage2;
Image1=double(OriginImage1)/256;
Image2=double(OriginImage2)/256;
[c1,s1]=wavedec2(Image1,2,'sym3'); %将X1进行2维分解，并使用sym4小波进行变换
[c2,s2]=wavedec2(Image2,2,'sym3');
c=0.5*(c1+c2); %计算系数平均值
s=0.5*(s1+s2);
X=waverec2(c,s,'sym3'); %进行小波重构
handles.X=X;
guidata(hObject,handles);
axes(handles.axes3);
imshow(X);

%评价指标
X=double(X); %对融合图像预处理
[C,R]=size(X);

%相对标准差
s=size(size(Image1));
if s(2)==3 %判断是灰度图像还是RGB彩色图像
f1=rgb2gray(Image1);
f2=rgb2gray(Image2);
else
    f1=Image1;
    f2=Image2;
end 
G1=double(f1);
G2=double(f2);
[m1,n1]=size(G1);
[m2,n2]=size(G2);
u1=(sum(G1(:)))/(m1*n1);
u2=(sum(G2(:)))/(m2*n2);
c1=0;
c2=0;
for i=1:m1
    for j=1:n1
        w1=G1(i,j)-u1;
        w2=G2(i,j)-u2;
        c1=c1+w1^2;
        c2=c2+w2^2;
    end
end
f1=sqrt(c1/(m1*n1));
f2=sqrt(c2/(m2*n2));
f=(f1-f2)/f1;
set(handles.edit1,'String',num2str(f));


%峰值信噪比
MAX=255; %图像灰度级最大值
D=Image1-X;
MES=sum(D(:).*D(:))./prod(size(Image1));
PSNR=10*log10(MAX^2/sqrt(MES));
set(handles.edit2,'String',num2str(PSNR));


%空间频率
RF=0;
CF=0;
for fi=1:C-1
    for fj=1:R-1
        RF=RF+(X(fi,fj)-X(fi,fj+1)).^2;
    end
end
RF=sqrt(RF/(C*R));
for fi=1:C-1
    for fj=1:R-1
        CF=CF+(X(fi,fj)-X(fi+1,fj)).^2;
    end
end
CF=sqrt(CF/(C*R));
SF=sqrt(RF+CF);
set(handles.edit3,'String',num2str(SF));


%图像清晰度
n=C*R;
m=1;
for i=1:(C-1)
    for j=1:(R-1)
        x=X(i,j)-X(i,j+1);
        y=X(i,j)-X(i+1,j);
        z(m,1)=sqrt((x.^2+y.^2)/2);
        m=m+1;
    end
end
G=sum(z)/n;                                  
set(handles.edit4,'String',num2str(G));


%互信息
s1=size(size(X));
if s1(2)==3 %判断是灰度图像还是RGB彩色图像
    a=rgb2gray(OriginImage1);
    a=double(a);
    b=rgb2gray(OriginImage2);
    b=double(b);
else
    a=double(OriginImage1);
    b=double(OriginImage2);
end
[Ma,Na] = size(a);
[Mb,Nb] = size(b);
M=min(Ma,Mb);
N=min(Na,Nb);

%初始化直方图数组
hab = zeros(256,256);
ha = zeros(1,256);
hb = zeros(1,256);

%归一化
if max(max(a))~=min(min(a))
    a = (a-min(min(a)))./(max(max(a))-min(min(a)));
else
    a = zeros(M,N);
end

if max(max(b))-min(min(b))
    b = (b-min(min(b)))./(max(max(b))-min(min(b)));
else
    b = zeros(M,N);
end

a = double(int16(a*255))+1;
b = double(int16(b*255))+1;

%统计直方图
for i=1:M
    for j=1:N
       indexx =  a(i,j);
       indexy = b(i,j) ;
       hab(indexx,indexy) = hab(indexx,indexy)+1; %联合直方图
       ha(indexx) = ha(indexx)+1; %a图直方图
       hb(indexy) = hb(indexy)+1; %b图直方图
   end
end

%计算联合信息熵
hsum = sum(sum(hab));
index = find(hab~=0);
p = hab/hsum;
Hab = sum(sum(-p(index).*log(p(index))));

%计算a图信息熵
hsum = sum(sum(ha));
index = find(ha~=0);
p = ha/hsum;
Ha = sum(sum(-p(index).*log(p(index))));

%计算b图信息熵
hsum = sum(sum(hb));
index = find(hb~=0);
p = hb/hsum;
Hb = sum(sum(-p(index).*log(p(index))));

%计算a图和b图的互信息
MI = Ha+Hb-Hab;
set(handles.edit5,'String',num2str(MI));


%交叉熵
s=size(size(X));
if s(2)==3 %判断是灰度图像还是RGB彩色图像
f1=rgb2gray(OriginImage1);
f2=rgb2gray(OriginImage2);
else
    f1=OriginImage1;
    f2=OriginImage2;
end 
G1=double(f1);
G2=double(f2);
[m1,n1]=size(G1);
[m2,n2]=size(G2);
m2=m1;
n2=n1;
X1=zeros(1,256);
X2=zeros(1,256);
result=0;

%统计两图各灰度级像素
for i=1:m1
    for j=1:n1
        X1(G1(i,j)+1)=X1(G1(i,j)+1)+1;
        X2(G2(i,j)+1)=X2(G2(i,j)+1)+1;
    end
end

%计算两图各灰度级像素出现的概率
for k=1:256
    P1(k)=X1(k)/(m1*n1);
    P2(k)=X2(k)/(m1*n1);
    if((P1(k)~=0)&&(P2(k)~=0))
        result=P1(k)*log2(P1(k)/P2(k))+result;
    end
end
f=result;
set(handles.edit6,'String',num2str(f));


function edit1_Callback(hObject, eventdata, handles)
% hObject    handle to edit1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit1 as text
%        str2double(get(hObject,'String')) returns contents of edit1 as a double


% --- Executes during object creation, after setting all properties.
function edit1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit2_Callback(hObject, eventdata, handles)
% hObject    handle to edit2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit2 as text
%        str2double(get(hObject,'String')) returns contents of edit2 as a double


% --- Executes during object creation, after setting all properties.
function edit2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit3_Callback(hObject, eventdata, handles)
% hObject    handle to edit3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit3 as text
%        str2double(get(hObject,'String')) returns contents of edit3 as a double


% --- Executes during object creation, after setting all properties.
function edit3_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit4_Callback(hObject, eventdata, handles)
% hObject    handle to edit4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit4 as text
%        str2double(get(hObject,'String')) returns contents of edit4 as a double


% --- Executes during object creation, after setting all properties.
function edit4_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit5_Callback(hObject, eventdata, handles)
% hObject    handle to edit5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit5 as text
%        str2double(get(hObject,'String')) returns contents of edit5 as a double


% --- Executes during object creation, after setting all properties.
function edit5_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit6_Callback(hObject, eventdata, handles)
% hObject    handle to edit2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit2 as text
%        str2double(get(hObject,'String')) returns contents of edit2 as a double


% --- Executes during object creation, after setting all properties.
function edit6_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
