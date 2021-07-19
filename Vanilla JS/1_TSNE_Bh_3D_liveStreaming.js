let DesiredPerplexity = 5;
let numberOfIterations = 1;
let LearningRatio = 1;
let Momentum = 0.8;
let TradeOff = 0.5;
let p = [];
let y = [];
let oldy = [];
let numberOfSamplesInX;
let Maximum;
let RangeX = 0;
let angle = 0;
let points = [];
const projection = [
    [1, 0, 0],
    [0, 1, 0],
];
function setup() {
    createCanvas(400, 400);
    LearningRatio = LearningRatio * 4; //By definition of dydt ... It doesn't make sense having it inside the loop.
    TradeOff = TradeOff * 1.73; //sqrt of two. Relationship between de side lenght and the diagonal of the square quadtree ... It doesn't make sense having it inside the loop.
    numberOfSamplesInX = X.length;
    let numberOfDimentions = X[1].length;
    y=zeros([numberOfSamplesInX-1,2]);
    oldy=zeros([numberOfSamplesInX-1,2]);
    //Compute Non Symetric Pair Wise Affinities [pj|i] with perplexity [Per]
    //pj|i=exp(-||xi-xj||^2/2Sigma^2)/Sum(exp(-||xi-xk||^2/2Sigma^2))
    //Perp(Pi)=2^H(Pi)
    //H(Pi)=-Sum(pj|i)*log2(pj|i)
    //pj|i and later will be pji
    let top1 = Array(numberOfSamplesInX).fill(0);
    let top2 = Array(numberOfSamplesInX).fill(0);
    let middle1 = Array(numberOfSamplesInX).fill(0);
    let middle2 = Array(numberOfSamplesInX).fill(0);
    let bottom1 = Array(numberOfSamplesInX).fill(0);
    let bottom2 = Array(numberOfSamplesInX).fill(0);
    p = []; //pj|i
    let aux;
    //
    let IndexElements=[];
    for(let i = 0; i < X.length; i++){
        IndexElements[i] = i;
    }
    let VantagePoint = new VantagePointElement();
    VantagePoint.SelectASeedAndFindMu(X, IndexElements);
    VantagePoint.SearchKNeighbors(X, numberOfSamplesInX, DesiredPerplexity);
    //
    for(let i = 0; i < numberOfSamplesInX; i++){
        top1[i] = -0.0001;
        bottom1[i] = -0.0001;
    }
    let IShouldStay;
    aux = getMeThePairWiseAffinities1(X, numberOfSamplesInX, numberOfDimentions);
    for(let iter = 0; iter <= 100; iter++){
        IShouldStay = true;
        p = getMeThePairWiseAffinities2(aux, numberOfSamplesInX, top1);
        top2 = GetMeThePerplexity(p, numberOfSamplesInX);
        for(let i = 0; i <= numberOfSamplesInX - 1; i++){
            if(top2[i] < DesiredPerplexity){
                IShouldStay = false;
                bottom1[i] = top1[i];
                top1[i] = (top1[i] - 1) * 2;
            }
        }
        if(IShouldStay){iter = 101};
    }
    for(let i = 0; i <= numberOfSamplesInX - 1; i++){
        middle1[i] = (top1[i] + bottom1[i]) / 2;
    }
    p = getMeThePairWiseAffinities2(aux, numberOfSamplesInX, middle1);
    middle2 = GetMeThePerplexity(p, numberOfSamplesInX);
    for(let iter = 0; iter <= 100; iter++){
        //Decision Maker (you can do better than this, see it later)
        for(let i = 0; i <= numberOfSamplesInX - 1; i++){
            if(Math.abs(middle2[i] - DesiredPerplexity) < 0.01){
            }else if(middle2[i] > DesiredPerplexity){
                top1[i] = middle1[i];
                top2[i] = middle2[i];
                middle1[i] = (top1[i] + bottom1[i]) / 2;
            }else{
                bottom1[i] = middle1[i];
                bottom2[i] = middle2[i];
                middle1[i] = (top1[i] + bottom1[i]) / 2;
            }
        }
        p = getMeThePairWiseAffinities2(aux, numberOfSamplesInX, middle1);
        middle2 = GetMeThePerplexity(p, numberOfSamplesInX);
    }
    //Set pij= ( pj|i + pi|j ) / 2 * n
    for(let i = 0; i <= numberOfSamplesInX - 2; i++){
        for(let z = 0; z <= VantagePointQueryArray[i].length - 1; z++){
            j = VantagePointQueryArray[i][z] - i - 1;
            if(j <= numberOfSamplesInX - i - 2 && j > 0){
                p[i][z] = p[i][z] / numberOfSamplesInX;
            }
        }
    }
    //Sample Initial Solution Y
    for(let i = 0; i <= numberOfSamplesInX - 2; i++){
        for(let n = 0; n <= 2; n++){
            y[i][n] = Math.random()-0.5;
        }
    }
    //for t=1 to T do
    //compute low-dimensional affinities qij
    //qij = (1+||yi-yj||^2)^-1) / Sum((1+||yk-yl||^2)^-1)
    //Compute gradient dCdY
    //dCdY=4*Sum((pij-qij)*(yi-yj)*(1+||yi-yj||^2))^-1
    //Set y(t)=y(t-1) + n dCdy + a(t) * (y(t-1)-y(t-2))
    //end
}
function draw(){
    background(0);
    YUpload(p, y, oldy, numberOfSamplesInX, numberOfIterations, Momentum, LearningRatio, VantagePointQueryArray, TradeOff);
    Maximum = 0;
    for(let i = 0; i < X.length; i++){
        if(Math.abs(y[i][0])>Maximum){Maximum = Math.abs(y[i][0])}
        if(Math.abs(y[i][1])>Maximum){Maximum = Math.abs(y[i][1])}
        if(Math.abs(y[i][2])>Maximum){Maximum = Math.abs(y[i][2])}
    }
    for(let i = 0; i < y.length; i++){
        points[i] = createVector(y[i][0]/Maximum, y[i][1]/Maximum, y[i][2]/Maximum);
    }
    translate(width / 2, height / 2);
    const rotationZ = [
        [cos(angle), -sin(angle), 0],
        [sin(angle), cos(angle), 0],
        [0, 0, 1],
    ];
    const rotationX = [
        [1, 0, 0],
        [0, cos(angle), -sin(angle)],
        [0, sin(angle), cos(angle)],
    ];
    const rotationY = [
        [cos(angle), 0, sin(angle)],
        [0, 1, 0],
        [-sin(angle), 0, cos(angle)],
    ];
    let projected = [];
    for (let i = 0; i < points.length; i++) {
        let rotated = matmul(rotationY, points[i]);
        rotated = matmul(rotationX, rotated);
        rotated = matmul(rotationZ, rotated);
        let projected2d = matmul(projection, rotated);
        projected2d.mult(200);
        projected[i] = projected2d;
        //point(projected2d.x, projected2d.y);
    }
    for (let i = 0; i < projected.length; i++) {
        strokeWeight(5);
        noFill();
        const v = projected[i];
        stroke(255);
        point(v.x, v.y);
    }
    angle += 0.005;
}
function getMeThePairWiseAffinities1(X, numberOfSamplesInX, numberOfDimentions){
    //This is the part of the code that is independant of minusTwoSigmaSquared.
    //I can play this part of the code only once during the sigma search.
    let p = Array(numberOfSamplesInX).fill(0); //pj|i
    for(let i = 0; i <= numberOfSamplesInX - 2; i++){
        p[i] = Array(VantagePointQueryArray[i].length).fill(0);
        for(let z = 0; z <= VantagePointQueryArray[i].length - 1; z++){
            j = VantagePointQueryArray[i][z] - i - 1;
            if(j <= numberOfSamplesInX - i - 2 && j > 0){
                for(let n = 0; n <= numberOfDimentions - 1; n++){
                    p[i][z] = p[i][z] + Math.pow(X[i][n] - X[j + i + 1][n], 2);
                }
            }
        }
    }
    return p;
}
function getMeThePairWiseAffinities2(auxiliar, numberOfSamplesInX, minusTwoSigmaSquared){
    // get the sum of pair wise affinities and the non normilized pair wise affinities
    let sumOfPairWiseAffinities = Array(numberOfSamplesInX).fill(0);
    let p = Array(numberOfSamplesInX).fill(0); //pj|i
    for(let i = 0; i <= numberOfSamplesInX - 2; i++){
        p[i] = Array(VantagePointQueryArray[i].length).fill(0);
        for(let z = 0; z <= VantagePointQueryArray[i].length - 1; z++){
            j = VantagePointQueryArray[i][z] - i - 1;
            if(j <= numberOfSamplesInX - i - 2 && j > 0){
                p[i][z] = Math.exp(auxiliar[i][z] / (minusTwoSigmaSquared[i] + 0.000001));
                sumOfPairWiseAffinities[i] = sumOfPairWiseAffinities[i] + p[i][z];
                sumOfPairWiseAffinities[j + i + 1] = sumOfPairWiseAffinities[j + i + 1] + p[i][z];
            }
        }
    }
    // get the normilized pair wise affinities
    for(let i = 0; i <= numberOfSamplesInX - 2; i++){
        for(let z = 0; z <= VantagePointQueryArray[i].length - 1; z++){
            j = VantagePointQueryArray[i][z] - i - 1;
            if(j <= numberOfSamplesInX - i - 2 && j > 0){
                p[i][z] = p[i][z] / (sumOfPairWiseAffinities[j + i + 1] + 0.000001);
            }
        }
    }
    return p;
}
function GetMeThePerplexity(p, numberOfSamplesInX){
    let Perplexities=Array(numberOfSamplesInX).fill(0);
    for(i = 0; i <= numberOfSamplesInX - 2; i++){
        for(let z = 0; z <= VantagePointQueryArray[i].length - 1; z++){
            j = VantagePointQueryArray[i][z] - i - 1;
            if(j <= numberOfSamplesInX - i - 2 && j > 0){
                let aux = p[i][z] * Math.log2(p[i][z] + 0.001);
                Perplexities[i] = Perplexities[i] + aux;
                Perplexities[j + i + 1] = Perplexities[j + i + 1] + aux;
            }
        }
    }
    for(i = 0; i <= numberOfSamplesInX - 1; i++){
        Perplexities[i] = Math.pow( 2, -Perplexities[i]);
    }
    return Perplexities;
}
function YUpload(p, y, oldy, numberOfSamplesInX, numberOfIterations, Momentum, LearningRatio, VantagePointQueryArray, TradeOff){
    let aux;
    let aux1;
    let j;
    let BiggestY = 2;
    let IndexElements=[];
    let ResultOT = new OctTreeResults();
    for(let i = 0; i < y.length; i++){
        IndexElements[i] = i;
    }
    for(let iter = 0; iter <= numberOfIterations; iter++){
        let Fattr = zeros([numberOfSamplesInX-1,2]);
        let Frep = zeros([numberOfSamplesInX-1,2]);
        let Sumq = 0;
        OctTree = new OctTreeElement([0,0,0], BiggestY);
        OctTree.InsertInBoxes(y, IndexElements);
        BiggestY = 0;
        //F attr
        for(let i = 0; i <= numberOfSamplesInX - 2; i++){
            for(let z = 0; z <= VantagePointQueryArray[i].length - 1; z++){
                j = VantagePointQueryArray[i][z] - i - 1;
                if(j <= numberOfSamplesInX - i - 2 && j > 0){
                    aux1=p[i][z] * CalculateZQij( i, j + i + 1);
                    for(let n = 0; n <= 2; n++){
                        aux = aux1 * (y[i][n] - y[j+ i + 1][n]);
                        Fattr[i][n] = Fattr[i][n] + aux;
                        Fattr[j + i + 1][n] = Fattr[j + i + 1][n] - aux;
                    }
                }
            }
        }
        for(let i = 0; i <= numberOfSamplesInX - 2; i++){
            OctTree.ListOfEquivalentBodiesOfI(y, i, TradeOff, ResultOT);
            for(let z = 0; z <= ResultOT.ResultOfTheQueryOT1.length - 1; z++){
                aux = CalculateZQij2( i, ResultOT.ResultOfTheQueryOT1[z]);
                Sumq = Sumq + ResultOT.ResultOfTheQueryOT2[z] * 2 / aux;
                aux1 = Math.pow(aux, 2);
                for(let n = 0; n <= 2; n++){
                    aux = aux1 * (y[i][n] - ResultOT.ResultOfTheQueryOT1[z][n] * ResultOT.ResultOfTheQueryOT2[z])
                    Frep[i][n] = Frep[i][n] - aux;
                }
            }
            for(let z = 0; z <= ResultOT.ResultOfTheQueryOT3.length - 1; z++){
                aux = CalculateZQij( i, ResultOT.ResultOfTheQueryOT3[z]);
                Sumq = Sumq +  2 / aux;
                aux1 = Math.pow(aux, 2);
                for(let n = 0; n <= 2; n++){
                    aux = aux1 * (y[i][n] - y[ResultOT.ResultOfTheQueryOT3[z]][n]);
                    Frep[i][n] = Frep[i][n] - aux;
                }
            }
        }
        //y adjustment
        for(let i = 0; i <= numberOfSamplesInX - 1; i++){
            for(let n = 0; n <= 2; n++){
                aux = y[i][n];
                y[i][n] = y[i][n] - LearningRatio * (Fattr[i][n] + Frep[i][n] / Sumq) + Momentum * (y[i][n] - oldy[i][n]);
                oldy[i][n] = aux;
                if(BiggestY < Math.abs(y[i][n])){BiggestY = Math.abs(y[i][n])}
            }
        }
    }
    function CalculateZQij( i, j){
        let aux=0;
        for(let n = 0; n <= 2; n++){
            aux = aux + Math.pow(y[i][n] - y[j][n], 2);
        }
        return Math.pow(1 + aux, -1);
    }
    function CalculateZQij2( i, CenterOfMass){
        let aux=0;
        for(let n = 0; n <= 2; n++){
            aux = aux + Math.pow(y[i][n] - CenterOfMass[n], 2);
        }
        return Math.pow(1 + aux, -1);
    }
}
function zeros(dimensions){
    var array = [];
    for (var i = 0; i <= dimensions[0]; ++i) {
        array.push(dimensions.length == 1 ? 0 : zeros(dimensions.slice(1)));
    }
    return array;
}
