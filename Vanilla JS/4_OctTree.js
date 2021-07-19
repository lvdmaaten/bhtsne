let OctTreeElementCapacity = 3;
let OctTreeMaximumInception = 8;
class OctTreeResults{
    constructor(){
        this.ResultOfTheQueryOT1 = [];
        this.ResultOfTheQueryOT2 = [];
        this.ResultOfTheQueryOT3 = [];
    }
}
class OctTreeElement{
    constructor(Center, Radius){
        this.Center = Center;
        this.Radius = Radius;
        this.ListOfIndexElementsULB = [];
        this.ListOfIndexElementsURB = [];
        this.ListOfIndexElementsDLB = [];
        this.ListOfIndexElementsDRB = [];
        this.ListOfIndexElementsULF = [];
        this.ListOfIndexElementsURF = [];
        this.ListOfIndexElementsDLF = [];
        this.ListOfIndexElementsDRF = [];
        this.OctTreeULBSon = undefined;
        this.OctTreeURBSon = undefined;
        this.OctTreeDLBSon = undefined;
        this.OctTreeDRBSon = undefined;
        this.OctTreeULFSon = undefined;
        this.OctTreeURFSon = undefined;
        this.OctTreeDLFSon = undefined;
        this.OctTreeDRFSon = undefined;
        this.AmountOfElementsInMe = undefined;
        this.CenterOfMass = undefined;
        this.DistanceToTheSample = undefined;
        this.level = 0;
    }
    ListOfEquivalentBodiesOfI(DataBase, i, TradeOff, Result){
        Result.ResultOfTheQueryOT1 = [];
        Result.ResultOfTheQueryOT2 = [];
        Result.ResultOfTheQueryOT3 = [];
        this.EquivalentBodies(DataBase, i, TradeOff, this.OctTreeULBSon, this.ListOfIndexElementsULB, Result);
        this.EquivalentBodies(DataBase, i, TradeOff, this.OctTreeURBSon, this.ListOfIndexElementsURB, Result);
        this.EquivalentBodies(DataBase, i, TradeOff, this.OctTreeDLBSon, this.ListOfIndexElementsDLB, Result);
        this.EquivalentBodies(DataBase, i, TradeOff, this.OctTreeDRBSon, this.ListOfIndexElementsDRB, Result);
        this.EquivalentBodies(DataBase, i, TradeOff, this.OctTreeULFSon, this.ListOfIndexElementsULF, Result);
        this.EquivalentBodies(DataBase, i, TradeOff, this.OctTreeURFSon, this.ListOfIndexElementsURF, Result);
        this.EquivalentBodies(DataBase, i, TradeOff, this.OctTreeDLFSon, this.ListOfIndexElementsDLF, Result);
        this.EquivalentBodies(DataBase, i, TradeOff, this.OctTreeDRFSon, this.ListOfIndexElementsDRF, Result);
    }
    EquivalentBodies(DataBase, i, TradeOff, OctTree, ListOfElements, Result){
        let aux;
        if(typeof OctTree != 'undefined'){
            aux = Math.pow(DataBase[i][0]-OctTree.CenterOfMass[0],2);
            aux = aux + Math.pow(DataBase[i][1]-OctTree.CenterOfMass[1],2);
            aux = aux + Math.pow(DataBase[i][2]-OctTree.CenterOfMass[2],2);
            aux = Math.pow(aux,0.5);
            if(OctTree.Radius / aux < TradeOff){
                Result.ResultOfTheQueryOT1.push(OctTree.CenterOfMass);
                Result.ResultOfTheQueryOT2.push(OctTree.AmountOfElementsInMe);
            }else{
                OctTree.EquivalentBodies(DataBase, i, TradeOff, OctTree.OctTreeULBSon, OctTree.ListOfIndexElementsULB, Result);
                OctTree.EquivalentBodies(DataBase, i, TradeOff, OctTree.OctTreeURBSon, OctTree.ListOfIndexElementsURB, Result);
                OctTree.EquivalentBodies(DataBase, i, TradeOff, OctTree.OctTreeDLBSon, OctTree.ListOfIndexElementsDLB, Result);
                OctTree.EquivalentBodies(DataBase, i, TradeOff, OctTree.OctTreeDRBSon, OctTree.ListOfIndexElementsDRB, Result);
                OctTree.EquivalentBodies(DataBase, i, TradeOff, OctTree.OctTreeULFSon, OctTree.ListOfIndexElementsULF, Result);
                OctTree.EquivalentBodies(DataBase, i, TradeOff, OctTree.OctTreeURFSon, OctTree.ListOfIndexElementsURF, Result);
                OctTree.EquivalentBodies(DataBase, i, TradeOff, OctTree.OctTreeDLFSon, OctTree.ListOfIndexElementsDLF, Result);
                OctTree.EquivalentBodies(DataBase, i, TradeOff, OctTree.OctTreeDRFSon, OctTree.ListOfIndexElementsDRF, Result);
            }
        }else{
            for(let i = 0; i < ListOfElements.length; i++){
                Result.ResultOfTheQueryOT3.push(ListOfElements[i]);
            }
        }
    }
    InsertInBoxes(DataBase, IndexOfElementsToEvaluate){
        this.AmountOfElementsInMe = IndexOfElementsToEvaluate.length;
        for(let Sample = 0; Sample < this.AmountOfElementsInMe; Sample++){
            if (DataBase[IndexOfElementsToEvaluate[Sample]][0] < this.Center[0]){
                if (DataBase[IndexOfElementsToEvaluate[Sample]][1] < this.Center[1]){
                    if (DataBase[IndexOfElementsToEvaluate[Sample]][2] < this.Center[2]){
                        this.ListOfIndexElementsDLB.push(IndexOfElementsToEvaluate[Sample]);
                    }else{
                        this.ListOfIndexElementsDLF.push(IndexOfElementsToEvaluate[Sample]);
                    }
                }else{
                    if (DataBase[IndexOfElementsToEvaluate[Sample]][2] < this.Center[2]){
                        this.ListOfIndexElementsULB.push(IndexOfElementsToEvaluate[Sample]);
                    }else{
                        this.ListOfIndexElementsULF.push(IndexOfElementsToEvaluate[Sample]);
                    }
                }
            }else{
                if (DataBase[IndexOfElementsToEvaluate[Sample]][1] < this.Center[1]){
                    if (DataBase[IndexOfElementsToEvaluate[Sample]][2] < this.Center[2]){
                        this.ListOfIndexElementsDRB.push(IndexOfElementsToEvaluate[Sample]);
                    }else{
                        this.ListOfIndexElementsDRF.push(IndexOfElementsToEvaluate[Sample]);
                    }
                }else{
                    if (DataBase[IndexOfElementsToEvaluate[Sample]][2] < this.Center[2]){
                        this.ListOfIndexElementsURB.push(IndexOfElementsToEvaluate[Sample]);
                    }else{
                        this.ListOfIndexElementsURF.push(IndexOfElementsToEvaluate[Sample]);
                    }
                }
            }
        }
        this.CalculateCenterOfMass(DataBase, IndexOfElementsToEvaluate);
        this.ShouldIBeDivided(DataBase);
    }
    CalculateCenterOfMass(DataBase, IndexOfElementsToEvaluate){
        let sumX = 0;
        let sumY = 0;
        let sumZ = 0;
        for(let Sample=0; Sample < this.AmountOfElementsInMe; Sample++){
            sumX = sumX + DataBase[IndexOfElementsToEvaluate[Sample]][0];
            sumY = sumY + DataBase[IndexOfElementsToEvaluate[Sample]][1];
            sumZ = sumZ + DataBase[IndexOfElementsToEvaluate[Sample]][2];
        }
        this.CenterOfMass = [sumX/this.AmountOfElementsInMe, sumY/this.AmountOfElementsInMe, sumZ/this.AmountOfElementsInMe];
    }
    ShouldIBeDivided(DataBase){
        if(this.ListOfIndexElementsULB.length > OctTreeElementCapacity && this.level < OctTreeMaximumInception){
            this.SubDivideULB(DataBase);
        }
        if(this.ListOfIndexElementsURB.length > OctTreeElementCapacity && this.level < OctTreeMaximumInception){
            this.SubDivideURB(DataBase);
        }
        if(this.ListOfIndexElementsDLB.length > OctTreeElementCapacity && this.level < OctTreeMaximumInception){
            this.SubDivideDLB(DataBase);
        }
        if(this.ListOfIndexElementsDRB.length > OctTreeElementCapacity && this.level < OctTreeMaximumInception){
            this.SubDivideDRB(DataBase);
        }
        if(this.ListOfIndexElementsULF.length > OctTreeElementCapacity && this.level < OctTreeMaximumInception){
            this.SubDivideULF(DataBase);
        }
        if(this.ListOfIndexElementsURF.length > OctTreeElementCapacity && this.level < OctTreeMaximumInception){
            this.SubDivideURF(DataBase);
        }
        if(this.ListOfIndexElementsDLF.length > OctTreeElementCapacity && this.level < OctTreeMaximumInception){
            this.SubDivideDLF(DataBase);
        }
        if(this.ListOfIndexElementsDRF.length > OctTreeElementCapacity && this.level < OctTreeMaximumInception){
            this.SubDivideDRF(DataBase);
        }
    }
    SubDivideULB(DataBase){
        let NewRadius = this.Radius/2;
        let NewCenter = [this.Center[0] - NewRadius,this.Center[1] + NewRadius,this.Center[2] - NewRadius];
        this.OctTreeULBSon = new OctTreeElement(NewCenter, NewRadius);
        this.OctTreeULBSon.level = this.level + 1;
        this.OctTreeULBSon.InsertInBoxes(DataBase, this.ListOfIndexElementsULB);
        //this.ListOfIndexElementsULB = []; //This is recommended to aliviate memory. but is not mandatory. Specially if you are looking for an error;
    }
    SubDivideURB(DataBase){
        let NewRadius = this.Radius/2;
        let NewCenter = [this.Center[0] + NewRadius,this.Center[1] + NewRadius,this.Center[2] - NewRadius];
        this.OctTreeURBSon = new OctTreeElement(NewCenter, NewRadius);
        this.OctTreeURBSon.level = this.level + 1;
        this.OctTreeURBSon.InsertInBoxes(DataBase, this.ListOfIndexElementsURB);
        //this.ListOfIndexElementsURB = []; //This is recommended to aliviate memory. but is not mandatory. Specially if you are looking for an error;
    }
    SubDivideDLB(DataBase){
        let NewRadius = this.Radius/2;
        let NewCenter = [this.Center[0] - NewRadius,this.Center[1] - NewRadius,this.Center[2] - NewRadius];
        this.OctTreeDLBSon = new OctTreeElement(NewCenter, NewRadius);
        this.OctTreeDLBSon.level = this.level + 1;
        this.OctTreeDLBSon.InsertInBoxes(DataBase, this.ListOfIndexElementsDLB);
        //this.ListOfIndexElementsDLB = []; //This is recommended to aliviate memory. but is not mandatory. Specially if you are looking for an error;
    }
    SubDivideDRB(DataBase){
        let NewRadius = this.Radius/2;
        let NewCenter = [this.Center[0] + NewRadius,this.Center[1] - NewRadius,this.Center[2] - NewRadius];
        this.OctTreeDRBSon = new OctTreeElement(NewCenter, NewRadius);
        this.OctTreeDRBSon.level = this.level + 1;
        this.OctTreeDRBSon.InsertInBoxes(DataBase, this.ListOfIndexElementsDRB);
        //this.ListOfIndexElementsDRB = []; //This is recommended to aliviate memory. but is not mandatory. Specially if you are looking for an error;
    }
    SubDivideULF(DataBase){
        let NewRadius = this.Radius/2;
        let NewCenter = [this.Center[0] - NewRadius,this.Center[1] + NewRadius,this.Center[2] + NewRadius];
        this.OctTreeULFSon = new OctTreeElement(NewCenter, NewRadius);
        this.OctTreeULFSon.level = this.level + 1;
        this.OctTreeULFSon.InsertInBoxes(DataBase, this.ListOfIndexElementsULF);
        //this.ListOfIndexElementsULF = []; //This is recommended to aliviate memory. but is not mandatory. Specially if you are looking for an error;
    }
    SubDivideURF(DataBase){
        let NewRadius = this.Radius/2;
        let NewCenter = [this.Center[0] + NewRadius,this.Center[1] + NewRadius,this.Center[2] + NewRadius];
        this.OctTreeURFSon = new OctTreeElement(NewCenter, NewRadius);
        this.OctTreeURFSon.level = this.level + 1;
        this.OctTreeURFSon.InsertInBoxes(DataBase, this.ListOfIndexElementsURF);
        //this.ListOfIndexElementsURF = []; //This is recommended to aliviate memory. but is not mandatory. Specially if you are looking for an error;
    }
    SubDivideDLF(DataBase){
        let NewRadius = this.Radius/2;
        let NewCenter = [this.Center[0] - NewRadius,this.Center[1] - NewRadius,this.Center[2] + NewRadius];
        this.OctTreeDLFSon = new OctTreeElement(NewCenter, NewRadius);
        this.OctTreeDLFSon.level = this.level + 1;
        this.OctTreeDLFSon.InsertInBoxes(DataBase, this.ListOfIndexElementsDLF);
        //this.ListOfIndexElementsDLF = []; //This is recommended to aliviate memory. but is not mandatory. Specially if you are looking for an error;
    }
    SubDivideDRF(DataBase){
        let NewRadius = this.Radius/2;
        let NewCenter = [this.Center[0] + NewRadius,this.Center[1] - NewRadius,this.Center[2] + NewRadius];
        this.OctTreeDRFSon = new OctTreeElement(NewCenter, NewRadius);
        this.OctTreeDRFSon.level = this.level + 1;
        this.OctTreeDRFSon.InsertInBoxes(DataBase, this.ListOfIndexElementsDRF);
        //this.ListOfIndexElementsDRF = []; //This is recommended to aliviate memory. but is not mandatory. Specially if you are looking for an error;
    }
}