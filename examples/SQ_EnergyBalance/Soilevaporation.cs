
using System;
using System.Collections.Generic;
using System.Linq;
using System.Xml;
using CRA.ModelLayer.MetadataTypes;
using CRA.ModelLayer.Core;
using CRA.ModelLayer.Strategy;
using System.Reflection;
using VarInfo=CRA.ModelLayer.Core.VarInfo;
using Preconditions=CRA.ModelLayer.Core.Preconditions;
using CRA.AgroManagement;       
using EnergyBalance.DomainClass;
namespace EnergyBalance.Strategies
{
    public class SoilEvaporation : IStrategyEnergyBalance
    {
        public SoilEvaporation()
        {
            ModellingOptions mo0_0 = new ModellingOptions();
            //Parameters
            List<VarInfo> _parameters0_0 = new List<VarInfo>();
            mo0_0.Parameters=_parameters0_0;

            //Inputs
            List<PropertyDescription> _inputs0_0 = new List<PropertyDescription>();
            PropertyDescription pd1 = new PropertyDescription();
            pd1.DomainClassType = typeof(EnergyBalance.DomainClass.EnergyBalanceState);
            pd1.PropertyName = "diffusionLimitedEvaporation";
            pd1.PropertyType = (EnergyBalance.DomainClass.EnergyBalanceStateVarInfo.diffusionLimitedEvaporation).ValueType.TypeForCurrentValue;
            pd1.PropertyVarInfo =(EnergyBalance.DomainClass.EnergyBalanceStateVarInfo.diffusionLimitedEvaporation);
            _inputs0_0.Add(pd1);
            PropertyDescription pd2 = new PropertyDescription();
            pd2.DomainClassType = typeof(EnergyBalance.DomainClass.EnergyBalanceState);
            pd2.PropertyName = "energyLimitedEvaporation";
            pd2.PropertyType = (EnergyBalance.DomainClass.EnergyBalanceStateVarInfo.energyLimitedEvaporation).ValueType.TypeForCurrentValue;
            pd2.PropertyVarInfo =(EnergyBalance.DomainClass.EnergyBalanceStateVarInfo.energyLimitedEvaporation);
            _inputs0_0.Add(pd2);
            mo0_0.Inputs=_inputs0_0;

            //Outputs
            List<PropertyDescription> _outputs0_0 = new List<PropertyDescription>();
            PropertyDescription pd3 = new PropertyDescription();
            pd3.DomainClassType = typeof(EnergyBalance.DomainClass.EnergyBalanceAuxiliary);
            pd3.PropertyName = "soilEvaporation";
            pd3.PropertyType = (EnergyBalance.DomainClass.EnergyBalanceAuxiliaryVarInfo.soilEvaporation).ValueType.TypeForCurrentValue;
            pd3.PropertyVarInfo =(EnergyBalance.DomainClass.EnergyBalanceAuxiliaryVarInfo.soilEvaporation);
            _outputs0_0.Add(pd3);
            mo0_0.Outputs=_outputs0_0;
            //Associated strategies
            List<string> lAssStrat0_0 = new List<string>();
            mo0_0.AssociatedStrategies = lAssStrat0_0;
            //Adding the modeling options to the modeling options manager
            _modellingOptionsManager = new ModellingOptionsManager(mo0_0);
            SetStaticParametersVarInfoDefinitions();
            SetPublisherData();

        }

        public string Description
        {
            get { return "Starting from a soil at field capacity, soil evaporation  is assumed to                be energy limited during the first phase of evaporation and diffusion limited thereafter.                Hence, the soil evaporation model considers these two processes taking the minimum between                the energy limited evaporation (PtSoil) and the diffused limited                evaporation " ;}
        }

        public string URL
        {
            get { return "" ;}
        }

        public string Domain
        {
            get { return "";}
        }

        public string ModelType
        {
            get { return "";}
        }

        public bool IsContext
        {
            get { return false;}
        }

        public IList<int> TimeStep
        {
            get
            {
                IList<int> ts = new List<int>();
                return ts;
            }
        }

        private  PublisherData _pd;
        public PublisherData PublisherData
        {
            get { return _pd;} 
        }

        private  void SetPublisherData()
        {
            _pd = new CRA.ModelLayer.MetadataTypes.PublisherData();
            _pd.Add("Creator", "Peter D. Jamieson, Glen S. Francis, Derick R. Wilson, Robert J. Martin");
            _pd.Add("Date", "");
            _pd.Add("Publisher", "New Zealand Institute for Crop and Food Research Ltd.,
            New Zealand Institute for Crop and Food Research Ltd.,
            New Zealand Institute for Crop and Food Research Ltd.,
            New Zealand Institute for Crop and Food Research Ltd.
        ");
        }

        private ModellingOptionsManager _modellingOptionsManager;
        public ModellingOptionsManager ModellingOptionsManager
        {
            get { return _modellingOptionsManager; } 
        }

        public IEnumerable<Type> GetStrategyDomainClassesTypes()
        {
            return new List<Type>() {  typeof(EnergyBalance.DomainClass.EnergyBalanceState),  typeof(EnergyBalance.DomainClass.EnergyBalanceState), typeof(EnergyBalance.DomainClass.EnergyBalanceRate), typeof(EnergyBalance.DomainClass.EnergyBalanceAuxiliary), typeof(EnergyBalance.DomainClass.EnergyBalanceExogenous)};
        }

        // Getter and setters for the value of the parameters of the strategy. The actual parameters are stored into the ModelingOptionsManager of the strategy.


        public void SetParametersDefaultValue()
        {
            _modellingOptionsManager.SetParametersDefaultValue();
        }

        private static void SetStaticParametersVarInfoDefinitions()
        {
        }

        public string TestPostConditions(EnergyBalance.DomainClass.EnergyBalanceState s,EnergyBalance.DomainClass.EnergyBalanceState s1,EnergyBalance.DomainClass.EnergyBalanceRate r,EnergyBalance.DomainClass.EnergyBalanceAuxiliary a,EnergyBalance.DomainClass.EnergyBalanceExogenous ex,string callID)
        {
            try
            {
                //Set current values of the outputs to the static VarInfo representing the output properties of the domain classes
                EnergyBalance.DomainClass.EnergyBalanceAuxiliaryVarInfo.soilEvaporation.CurrentValue=a.soilEvaporation;
                ConditionsCollection prc = new ConditionsCollection();
                Preconditions pre = new Preconditions(); 
                RangeBasedCondition r3 = new RangeBasedCondition(EnergyBalance.DomainClass.EnergyBalanceAuxiliaryVarInfo.soilEvaporation);
                if(r3.ApplicableVarInfoValueTypes.Contains( EnergyBalance.DomainClass.EnergyBalanceAuxiliaryVarInfo.soilEvaporation.ValueType)){prc.AddCondition(r3);}
                string postConditionsResult = pre.VerifyPostconditions(prc, callID); if (!string.IsNullOrEmpty(postConditionsResult)) { pre.TestsOut(postConditionsResult, true, "PostConditions errors in strategy " + this.GetType().Name); } return postConditionsResult;
            }
            catch (Exception exception)
            {
                string msg = ".EnergyBalance, " + this.GetType().Name + ": Unhandled exception running post-condition test. ";
                throw new Exception(msg, exception);
            }
        }

        public string TestPreConditions(EnergyBalance.DomainClass.EnergyBalanceState s,EnergyBalance.DomainClass.EnergyBalanceState s1,EnergyBalance.DomainClass.EnergyBalanceRate r,EnergyBalance.DomainClass.EnergyBalanceAuxiliary a,EnergyBalance.DomainClass.EnergyBalanceExogenous ex,string callID)
        {
            try
            {
                //Set current values of the inputs to the static VarInfo representing the inputs properties of the domain classes
                EnergyBalance.DomainClass.EnergyBalanceStateVarInfo.diffusionLimitedEvaporation.CurrentValue=s.diffusionLimitedEvaporation;
                EnergyBalance.DomainClass.EnergyBalanceStateVarInfo.energyLimitedEvaporation.CurrentValue=s.energyLimitedEvaporation;
                ConditionsCollection prc = new ConditionsCollection();
                Preconditions pre = new Preconditions(); 
                RangeBasedCondition r1 = new RangeBasedCondition(EnergyBalance.DomainClass.EnergyBalanceStateVarInfo.diffusionLimitedEvaporation);
                if(r1.ApplicableVarInfoValueTypes.Contains( EnergyBalance.DomainClass.EnergyBalanceStateVarInfo.diffusionLimitedEvaporation.ValueType)){prc.AddCondition(r1);}
                RangeBasedCondition r2 = new RangeBasedCondition(EnergyBalance.DomainClass.EnergyBalanceStateVarInfo.energyLimitedEvaporation);
                if(r2.ApplicableVarInfoValueTypes.Contains( EnergyBalance.DomainClass.EnergyBalanceStateVarInfo.energyLimitedEvaporation.ValueType)){prc.AddCondition(r2);}
                string preConditionsResult = pre.VerifyPreconditions(prc, callID); if (!string.IsNullOrEmpty(preConditionsResult)) { pre.TestsOut(preConditionsResult, true, "PreConditions errors in strategy " + this.GetType().Name); } return preConditionsResult;
            }
            catch (Exception exception)
            {
                string msg = ".EnergyBalance, " + this.GetType().Name + ": Unhandled exception running pre-condition test. ";
                throw new Exception(msg, exception);
            }
        }

        public void Estimate(EnergyBalance.DomainClass.EnergyBalanceState s,EnergyBalance.DomainClass.EnergyBalanceState s1,EnergyBalance.DomainClass.EnergyBalanceRate r,EnergyBalance.DomainClass.EnergyBalanceAuxiliary a,EnergyBalance.DomainClass.EnergyBalanceExogenous ex)
        {
            try
            {
                CalculateModel(s, s1, r, a, ex);
            }
            catch (Exception exception)
            {
                string msg = "Error in component EnergyBalance, strategy: " + this.GetType().Name + ": Unhandled exception running model. "+exception.GetType().FullName+" - "+exception.Message;
                throw new Exception(msg, exception);
            }
        }

        private void CalculateModel(EnergyBalance.DomainClass.EnergyBalanceState s, EnergyBalance.DomainClass.EnergyBalanceState s1, EnergyBalance.DomainClass.EnergyBalanceRate r, EnergyBalance.DomainClass.EnergyBalanceAuxiliary a, EnergyBalance.DomainClass.EnergyBalanceExogenous ex)
        {
            double diffusionLimitedEvaporation = s.diffusionLimitedEvaporation;
            double energyLimitedEvaporation = s.energyLimitedEvaporation;
            double soilEvaporation;
            soilEvaporation = Math.Min(diffusionLimitedEvaporation, energyLimitedEvaporation);
            a.soilEvaporation= soilEvaporation;
        }

    }
}