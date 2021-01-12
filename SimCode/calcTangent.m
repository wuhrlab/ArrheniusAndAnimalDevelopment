function tang = calcTangent(x, y, temp)

    dy = diff(y)./diff(x);
    tang = (x-x(temp))*dy(temp)+y(temp);

end