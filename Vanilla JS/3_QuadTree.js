let QuadTreeElementCapacity = 3;
let QuadTreeMaximumInception = 8;
class QuadTreeResults{
    constructor(){
        this.ResultOfTheQueryOT1 = [];
        this.ResultOfTheQueryOT2 = [];
        this.ResultOfTheQueryOT3 = [];
    }
}
class QuadtreeElement{
    constructor(Center, Radius){
        this.Center = Center;
        this.Radius = Radius;
        this.ListOfIndexElementsUL = [];
        this.ListOfIndexElementsUR = [];
        this.ListOfIndexElementsDL = [];
        this.ListOfIndexElementsDR = [];
        this.QuadtreeULSon = undefined;
        this.QuadtreeURSon = undefined;
        this.QuadtreeDLSon = undefined;
        this.QuadtreeDRSon = undefined;
        this.AmountOfElementsInMe = undefined;
        this.CenterOfMass = undefined;
        this.DistanceToTheSample = undefined;
        this.level = 0;
    }
    ListOfEquivalentBodiesOfI(DataBase, i, TradeOff, Result){
        Result.ResultOfTheQueryQT1 = [];
        Result.ResultOfTheQueryQT2 = [];
        Result.ResultOfTheQueryQT3 = [];
        this.EquivalentBodies(DataBase, i, TradeOff, this.QuadtreeULSon, this.ListOfIndexElementsUL, Result);
        this.EquivalentBodies(DataBase, i, TradeOff, this.QuadtreeURSon, this.ListOfIndexElementsUR, Result);
        this.EquivalentBodies(DataBase, i, TradeOff, this.QuadtreeDLSon, this.ListOfIndexElementsDL, Result);
        this.EquivalentBodies(DataBase, i, TradeOff, this.QuadtreeDRSon, this.ListOfIndexElementsDR, Result);
    }
    EquivalentBodies(DataBase, i, TradeOff, QuadTree, ListOfElements, Result){
        let aux;
        if(typeof QuadTree != 'undefined'){
            aux = Math.pow(DataBase[i][0]-QuadTree.CenterOfMass[0],2);
            aux = aux + Math.pow(DataBase[i][1]-QuadTree.CenterOfMass[1],2);
            aux = Math.pow(aux,0.5);
            if(QuadTree.Radius / aux < TradeOff){
                Result.ResultOfTheQueryQT1.push(QuadTree.CenterOfMass);
                Result.ResultOfTheQueryQT2.push(QuadTree.AmountOfElementsInMe);
            }else{
                QuadTree.EquivalentBodies(DataBase, i, TradeOff, QuadTree.QuadtreeULSon, QuadTree.ListOfIndexElementsUL, Result);
                QuadTree.EquivalentBodies(DataBase, i, TradeOff, QuadTree.QuadtreeURSon, QuadTree.ListOfIndexElementsUR, Result);
                QuadTree.EquivalentBodies(DataBase, i, TradeOff, QuadTree.QuadtreeDLSon, QuadTree.ListOfIndexElementsDL, Result);
                QuadTree.EquivalentBodies(DataBase, i, TradeOff, QuadTree.QuadtreeDRSon, QuadTree.ListOfIndexElementsDR, Result);
            }
        }else{
            for(let i = 0; i < ListOfElements.length; i++){
                Result.ResultOfTheQueryQT3.push(ListOfElements[i]);
            }
        }
    }
    InsertInBoxes(DataBase, IndexOfElementsToEvaluate){
        this.AmountOfElementsInMe = IndexOfElementsToEvaluate.length;
        for(let Sample = 0; Sample < this.AmountOfElementsInMe; Sample++){
            if (DataBase[IndexOfElementsToEvaluate[Sample]][0] < this.Center[0]){
                if (DataBase[IndexOfElementsToEvaluate[Sample]][1] < this.Center[1]){
                    this.ListOfIndexElementsDL.push(IndexOfElementsToEvaluate[Sample]);
                }else{
                    this.ListOfIndexElementsUL.push(IndexOfElementsToEvaluate[Sample]);
                }
            }else{
                if (DataBase[IndexOfElementsToEvaluate[Sample]][1] < this.Center[1]){
                    this.ListOfIndexElementsDR.push(IndexOfElementsToEvaluate[Sample]);
                }else{
                    this.ListOfIndexElementsUR.push(IndexOfElementsToEvaluate[Sample]);
                }
            }
        }
        this.CalculateCenterOfMass(DataBase, IndexOfElementsToEvaluate);
        this.ShouldIBeDivided(DataBase);
    }
    CalculateCenterOfMass(DataBase, IndexOfElementsToEvaluate){
        let sumX = 0;
        let sumY = 0;
        for(let Sample=0; Sample < this.AmountOfElementsInMe; Sample++){
            sumX = sumX + DataBase[IndexOfElementsToEvaluate[Sample]][0];
            sumY = sumY + DataBase[IndexOfElementsToEvaluate[Sample]][1];
        }
        this.CenterOfMass = [sumX/this.AmountOfElementsInMe, sumY/this.AmountOfElementsInMe];
    }
    ShouldIBeDivided(DataBase){
        if(this.ListOfIndexElementsUL.length > QuadTreeElementCapacity && this.level < QuadTreeMaximumInception){
            this.SubDivideUL(DataBase);
        }
        if(this.ListOfIndexElementsUR.length > QuadTreeElementCapacity && this.level < QuadTreeMaximumInception){
            this.SubDivideUR(DataBase);
        }
        if(this.ListOfIndexElementsDL.length > QuadTreeElementCapacity && this.level < QuadTreeMaximumInception){
            this.SubDivideDL(DataBase);
        }
        if(this.ListOfIndexElementsDR.length > QuadTreeElementCapacity && this.level < QuadTreeMaximumInception){
            this.SubDivideDR(DataBase);
        }
    }
    SubDivideUL(DataBase){
        let NewRadius = this.Radius/2;
        let NewCenter = [this.Center[0] - NewRadius,this.Center[1] + NewRadius];
        this.QuadtreeULSon = new QuadtreeElement(NewCenter, NewRadius);
        this.QuadtreeULSon.level = this.level + 1;
        this.QuadtreeULSon.InsertInBoxes(DataBase, this.ListOfIndexElementsUL);
        this.ListOfIndexElementsUL = []; //This is recommended to aliviate memory. but is not mandatory. Specially if you are looking for an error;
    }
    SubDivideUR(DataBase){
        let NewRadius = this.Radius/2;
        let NewCenter = [this.Center[0] + NewRadius,this.Center[1] + NewRadius];
        this.QuadtreeURSon = new QuadtreeElement(NewCenter, NewRadius);
        this.QuadtreeURSon.level = this.level + 1;
        this.QuadtreeURSon.InsertInBoxes(DataBase, this.ListOfIndexElementsUR);
        this.ListOfIndexElementsUR = []; //This is recommended to aliviate memory. but is not mandatory. Specially if you are looking for an error;
    }
    SubDivideDL(DataBase){
        let NewRadius = this.Radius/2;
        let NewCenter = [this.Center[0] - NewRadius,this.Center[1] - NewRadius];
        this.QuadtreeDLSon = new QuadtreeElement(NewCenter, NewRadius);
        this.QuadtreeDLSon.level = this.level + 1;
        this.QuadtreeDLSon.InsertInBoxes(DataBase, this.ListOfIndexElementsDL);
        this.ListOfIndexElementsDL = []; //This is recommended to aliviate memory. but is not mandatory. Specially if you are looking for an error;
    }
    SubDivideDR(DataBase){
        let NewRadius = this.Radius/2;
        let NewCenter = [this.Center[0] + NewRadius,this.Center[1] - NewRadius];
        this.QuadtreeDRSon = new QuadtreeElement(NewCenter, NewRadius);
        this.QuadtreeDRSon.level = this.level + 1;
        this.QuadtreeDRSon.InsertInBoxes(DataBase, this.ListOfIndexElementsDR);
        this.ListOfIndexElementsDR = []; //This is recommended to aliviate memory. but is not mandatory. Specially if you are looking for an error;
    }
}