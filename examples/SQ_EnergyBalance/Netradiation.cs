
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
    public class NetRadiation : IStrategyEnergyBalance
    {
        public NetRadiation()
        {
            ModellingOptions mo0_0 = new ModellingOptions();
            //Parameters
            List<VarInfo> _parameters0_0 = new List<VarInfo>();
            VarInfo v1 = new VarInfo();
            v1.DefaultValue = 0.23;
            v1.Description = "albedo Coefficient";
            v1.Id = 0;
            v1.MaxValue = 1;
            v1.MinValue = 0;
            v1.Name = "albedoCoefficient";
            v1.Size = 1;
            v1.Units = "";
            v1.URL = "";
            v1.VarType = CRA.ModelLayer.Core.VarInfo.Type.PARAMETER;
            v1.ValueType = VarInfoValueTypes.GetInstanceForName("Double");
            _parameters0_0.Add(v1);
            VarInfo v2 = new VarInfo();
            v2.DefaultValue = 4.903E-09;
            v2.Description = "stefan Boltzman constant";
            v2.Id = 0;
            v2.MaxValue = 1;
            v2.MinValue = 0;
            v2.Name = "stefanBoltzman";
            v2.Size = 1;
            v2.Units = "";
            v2.URL = "";
            v2.VarType = CRA.ModelLayer.Core.VarInfo.Type.PARAMETER;
            v2.ValueType = VarInfoValueTypes.GetInstanceForName("Double");
            _parameters0_0.Add(v2);
            VarInfo v3 = new VarInfo();
            v3.DefaultValue = 0;
            v3.Description = "elevation";
            v3.Id = 0;
            v3.MaxValue = 10000;
            v3.MinValue = -500;
            v3.Name = "elevation";
            v3.Size = 1;
            v3.Units = "m";
            v3.URL = "";
            v3.VarType = CRA.ModelLayer.Core.VarInfo.Type.PARAMETER;
            v3.ValueType = VarInfoValueTypes.GetInstanceForName("Double");
            _parameters0_0.Add(v3);
            mo0_0.Parameters=_parameters0_0;

            //Inputs
            List<PropertyDescription> _inputs0_0 = new List<PropertyDescription>();
            PropertyDescription pd1 = new PropertyDescription();
            pd1.DomainClassType = typeof(EnergyBalance.DomainClass.EnergyBalanceAuxiliary);
            pd1.PropertyName = "minTair";
            pd1.PropertyType = (EnergyBalance.DomainClass.EnergyBalanceAuxiliaryVarInfo.minTair).ValueType.TypeForCurrentValue;
            pd1.PropertyVarInfo =(EnergyBalance.DomainClass.EnergyBalanceAuxiliaryVarInfo.minTair);
            _inputs0_0.Add(pd1);
            PropertyDescription pd2 = new PropertyDescription();
            pd2.DomainClassType = typeof(EnergyBalance.DomainClass.EnergyBalanceAuxiliary);
            pd2.PropertyName = "maxTair";
            pd2.PropertyType = (EnergyBalance.DomainClass.EnergyBalanceAuxiliaryVarInfo.maxTair).ValueType.TypeForCurrentValue;
            pd2.PropertyVarInfo =(EnergyBalance.DomainClass.EnergyBalanceAuxiliaryVarInfo.maxTair);
            _inputs0_0.Add(pd2);
            PropertyDescription pd3 = new PropertyDescription();
            pd3.DomainClassType = typeof(EnergyBalance.DomainClass.EnergyBalanceAuxiliary);
            pd3.PropertyName = "solarRadiation";
            pd3.PropertyType = (EnergyBalance.DomainClass.EnergyBalanceAuxiliaryVarInfo.solarRadiation).ValueType.TypeForCurrentValue;
            pd3.PropertyVarInfo =(EnergyBalance.DomainClass.EnergyBalanceAuxiliaryVarInfo.solarRadiation);
            _inputs0_0.Add(pd3);
            PropertyDescription pd4 = new PropertyDescription();
            pd4.DomainClassType = typeof(EnergyBalance.DomainClass.EnergyBalanceAuxiliary);
            pd4.PropertyName = "vaporPressure";
            pd4.PropertyType = (EnergyBalance.DomainClass.EnergyBalanceAuxiliaryVarInfo.vaporPressure).ValueType.TypeForCurrentValue;
            pd4.PropertyVarInfo =(EnergyBalance.DomainClass.EnergyBalanceAuxiliaryVarInfo.vaporPressure);
            _inputs0_0.Add(pd4);
            PropertyDescription pd5 = new PropertyDescription();
            pd5.DomainClassType = typeof(EnergyBalance.DomainClass.EnergyBalanceAuxiliary);
            pd5.PropertyName = "extraSolarRadiation";
            pd5.PropertyType = (EnergyBalance.DomainClass.EnergyBalanceAuxiliaryVarInfo.extraSolarRadiation).ValueType.TypeForCurrentValue;
            pd5.PropertyVarInfo =(EnergyBalance.DomainClass.EnergyBalanceAuxiliaryVarInfo.extraSolarRadiation);
            _inputs0_0.Add(pd5);
            mo0_0.Inputs=_inputs0_0;

            //Outputs
            List<PropertyDescription> _outputs0_0 = new List<PropertyDescription>();
            PropertyDescription pd6 = new PropertyDescription();
            pd6.DomainClassType = typeof(EnergyBalance.DomainClass.EnergyBalanceAuxiliary);
            pd6.PropertyName = "netRadiation";
            pd6.PropertyType = (EnergyBalance.DomainClass.EnergyBalanceAuxiliaryVarInfo.netRadiation).ValueType.TypeForCurrentValue;
            pd6.PropertyVarInfo =(EnergyBalance.DomainClass.EnergyBalanceAuxiliaryVarInfo.netRadiation);
            _outputs0_0.Add(pd6);
            mo0_0.Outputs=_outputs0_0;PropertyDescription pd7 = new PropertyDescription();
            pd7.DomainClassType = typeof(EnergyBalance.DomainClass.EnergyBalanceAuxiliary);
            pd7.PropertyName = "netOutGoingLongWaveRadiation";
            pd7.PropertyType = (EnergyBalance.DomainClass.EnergyBalanceAuxiliaryVarInfo.netOutGoingLongWaveRadiation).ValueType.TypeForCurrentValue;
            pd7.PropertyVarInfo =(EnergyBalance.DomainClass.EnergyBalanceAuxiliaryVarInfo.netOutGoingLongWaveRadiation);
            _outputs0_0.Add(pd7);
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
            get { return "It is calculated at the surface of the canopy and is givenby the difference between incoming and outgoing radiation of both short                     and long wavelength radiation " ;}
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

        public double albedoCoefficient
        {
            get { 
                VarInfo vi= _modellingOptionsManager.GetParameterByName("albedoCoefficient");
                if (vi != null && vi.CurrentValue!=null) return (double)vi.CurrentValue ;
                else throw new Exception("Parameter 'albedoCoefficient' not found (or found null) in strategy 'NetRadiation'");
            } set {
                VarInfo vi = _modellingOptionsManager.GetParameterByName("albedoCoefficient");
                if (vi != null)  vi.CurrentValue=value;
                else throw new Exception("Parameter 'albedoCoefficient' not found in strategy 'NetRadiation'");
            }
        }
        public double stefanBoltzman
        {
            get { 
                VarInfo vi= _modellingOptionsManager.GetParameterByName("stefanBoltzman");
                if (vi != null && vi.CurrentValue!=null) return (double)vi.CurrentValue ;
                else throw new Exception("Parameter 'stefanBoltzman' not found (or found null) in strategy 'NetRadiation'");
            } set {
                VarInfo vi = _modellingOptionsManager.GetParameterByName("stefanBoltzman");
                if (vi != null)  vi.CurrentValue=value;
                else throw new Exception("Parameter 'stefanBoltzman' not found in strategy 'NetRadiation'");
            }
        }
        public double elevation
        {
            get { 
                VarInfo vi= _modellingOptionsManager.GetParameterByName("elevation");
                if (vi != null && vi.CurrentValue!=null) return (double)vi.CurrentValue ;
                else throw new Exception("Parameter 'elevation' not found (or found null) in strategy 'NetRadiation'");
            } set {
                VarInfo vi = _modellingOptionsManager.GetParameterByName("elevation");
                if (vi != null)  vi.CurrentValue=value;
                else throw new Exception("Parameter 'elevation' not found in strategy 'NetRadiation'");
            }
        }

        public void SetParametersDefaultValue()
        {
            _modellingOptionsManager.SetParametersDefaultValue();
        }

        private static void SetStaticParametersVarInfoDefinitions()
        {

            albedoCoefficientVarInfo.Name = "albedoCoefficient";
            albedoCoefficientVarInfo.Description = "albedo Coefficient";
            albedoCoefficientVarInfo.MaxValue = 1;
            albedoCoefficientVarInfo.MinValue = 0;
            albedoCoefficientVarInfo.DefaultValue = 0.23;
            albedoCoefficientVarInfo.Units = "";
            albedoCoefficientVarInfo.ValueType = VarInfoValueTypes.GetInstanceForName("Double");

            stefanBoltzmanVarInfo.Name = "stefanBoltzman";
            stefanBoltzmanVarInfo.Description = "stefan Boltzman constant";
            stefanBoltzmanVarInfo.MaxValue = 1;
            stefanBoltzmanVarInfo.MinValue = 0;
            stefanBoltzmanVarInfo.DefaultValue = 4.903E-09;
            stefanBoltzmanVarInfo.Units = "";
            stefanBoltzmanVarInfo.ValueType = VarInfoValueTypes.GetInstanceForName("Double");

            elevationVarInfo.Name = "elevation";
            elevationVarInfo.Description = "elevation";
            elevationVarInfo.MaxValue = 10000;
            elevationVarInfo.MinValue = -500;
            elevationVarInfo.DefaultValue = 0;
            elevationVarInfo.Units = "m";
            elevationVarInfo.ValueType = VarInfoValueTypes.GetInstanceForName("Double");
        }

        private static VarInfo _albedoCoefficientVarInfo = new VarInfo();
        public static VarInfo albedoCoefficientVarInfo
        {
            get { return _albedoCoefficientVarInfo;} 
        }

        private static VarInfo _stefanBoltzmanVarInfo = new VarInfo();
        public static VarInfo stefanBoltzmanVarInfo
        {
            get { return _stefanBoltzmanVarInfo;} 
        }

        private static VarInfo _elevationVarInfo = new VarInfo();
        public static VarInfo elevationVarInfo
        {
            get { return _elevationVarInfo;} 
        }

        public string TestPostConditions(EnergyBalance.DomainClass.EnergyBalanceState s,EnergyBalance.DomainClass.EnergyBalanceState s1,EnergyBalance.DomainClass.EnergyBalanceRate r,EnergyBalance.DomainClass.EnergyBalanceAuxiliary a,EnergyBalance.DomainClass.EnergyBalanceExogenous ex,string callID)
        {
            try
            {
                //Set current values of the outputs to the static VarInfo representing the output properties of the domain classes
                EnergyBalance.DomainClass.EnergyBalanceAuxiliaryVarInfo.netRadiation.CurrentValue=a.netRadiation;
                EnergyBalance.DomainClass.EnergyBalanceAuxiliaryVarInfo.netOutGoingLongWaveRadiation.CurrentValue=a.netOutGoingLongWaveRadiation;
                ConditionsCollection prc = new ConditionsCollection();
                Preconditions pre = new Preconditions(); 
                RangeBasedCondition r9 = new RangeBasedCondition(EnergyBalance.DomainClass.EnergyBalanceAuxiliaryVarInfo.netRadiation);
                if(r9.ApplicableVarInfoValueTypes.Contains( EnergyBalance.DomainClass.EnergyBalanceAuxiliaryVarInfo.netRadiation.ValueType)){prc.AddCondition(r9);}
                RangeBasedCondition r10 = new RangeBasedCondition(EnergyBalance.DomainClass.EnergyBalanceAuxiliaryVarInfo.netOutGoingLongWaveRadiation);
                if(r10.ApplicableVarInfoValueTypes.Contains( EnergyBalance.DomainClass.EnergyBalanceAuxiliaryVarInfo.netOutGoingLongWaveRadiation.ValueType)){prc.AddCondition(r10);}
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
                EnergyBalance.DomainClass.EnergyBalanceAuxiliaryVarInfo.minTair.CurrentValue=a.minTair;
                EnergyBalance.DomainClass.EnergyBalanceAuxiliaryVarInfo.maxTair.CurrentValue=a.maxTair;
                EnergyBalance.DomainClass.EnergyBalanceAuxiliaryVarInfo.solarRadiation.CurrentValue=a.solarRadiation;
                EnergyBalance.DomainClass.EnergyBalanceAuxiliaryVarInfo.vaporPressure.CurrentValue=a.vaporPressure;
                EnergyBalance.DomainClass.EnergyBalanceAuxiliaryVarInfo.extraSolarRadiation.CurrentValue=a.extraSolarRadiation;
                ConditionsCollection prc = new ConditionsCollection();
                Preconditions pre = new Preconditions(); 
                RangeBasedCondition r1 = new RangeBasedCondition(EnergyBalance.DomainClass.EnergyBalanceAuxiliaryVarInfo.minTair);
                if(r1.ApplicableVarInfoValueTypes.Contains( EnergyBalance.DomainClass.EnergyBalanceAuxiliaryVarInfo.minTair.ValueType)){prc.AddCondition(r1);}
                RangeBasedCondition r2 = new RangeBasedCondition(EnergyBalance.DomainClass.EnergyBalanceAuxiliaryVarInfo.maxTair);
                if(r2.ApplicableVarInfoValueTypes.Contains( EnergyBalance.DomainClass.EnergyBalanceAuxiliaryVarInfo.maxTair.ValueType)){prc.AddCondition(r2);}
                RangeBasedCondition r3 = new RangeBasedCondition(EnergyBalance.DomainClass.EnergyBalanceAuxiliaryVarInfo.solarRadiation);
                if(r3.ApplicableVarInfoValueTypes.Contains( EnergyBalance.DomainClass.EnergyBalanceAuxiliaryVarInfo.solarRadiation.ValueType)){prc.AddCondition(r3);}
                RangeBasedCondition r4 = new RangeBasedCondition(EnergyBalance.DomainClass.EnergyBalanceAuxiliaryVarInfo.vaporPressure);
                if(r4.ApplicableVarInfoValueTypes.Contains( EnergyBalance.DomainClass.EnergyBalanceAuxiliaryVarInfo.vaporPressure.ValueType)){prc.AddCondition(r4);}
                RangeBasedCondition r5 = new RangeBasedCondition(EnergyBalance.DomainClass.EnergyBalanceAuxiliaryVarInfo.extraSolarRadiation);
                if(r5.ApplicableVarInfoValueTypes.Contains( EnergyBalance.DomainClass.EnergyBalanceAuxiliaryVarInfo.extraSolarRadiation.ValueType)){prc.AddCondition(r5);}
                prc.AddCondition(new RangeBasedCondition(_modellingOptionsManager.GetParameterByName("albedoCoefficient")));
                prc.AddCondition(new RangeBasedCondition(_modellingOptionsManager.GetParameterByName("stefanBoltzman")));
                prc.AddCondition(new RangeBasedCondition(_modellingOptionsManager.GetParameterByName("elevation")));
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
            double minTair = a.minTair;
            double maxTair = a.maxTair;
            double solarRadiation = a.solarRadiation;
            double vaporPressure = a.vaporPressure;
            double extraSolarRadiation = a.extraSolarRadiation;
            double netRadiation;
            double netOutGoingLongWaveRadiation;
            double Nsr;
            double clearSkySolarRadiation;
            double averageT;
            double surfaceEmissivity;
            double cloudCoverFactor;
            double Nolr;
            Nsr = (1.0d - albedoCoefficient) * solarRadiation;
            clearSkySolarRadiation = (0.75d + (2 * Math.Pow(10.0d, -5) * elevation)) * extraSolarRadiation;
            averageT = (Math.Pow(maxTair + 273.16d, 4) + Math.Pow(minTair + 273.16d, 4)) / 2.0d;
            surfaceEmissivity = 0.34d - (0.14d * Math.Sqrt(vaporPressure / 10.0d));
            cloudCoverFactor = 1.35d * (solarRadiation / clearSkySolarRadiation) - 0.35d;
            Nolr = stefanBoltzman * averageT * surfaceEmissivity * cloudCoverFactor;
            netRadiation = Nsr - Nolr;
            netOutGoingLongWaveRadiation = Nolr;
            a.netRadiation= netRadiation;
            a.netOutGoingLongWaveRadiation= netOutGoingLongWaveRadiation;
        }

    }
}