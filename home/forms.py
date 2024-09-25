from django import forms

class FastaUploadForm(forms.Form):
    score_model_conservation = forms.ChoiceField(
        choices=[(1,'Aminoacid classification'), (2,'BLOSUM62')],
        widget=forms.RadioSelect(attrs={
            'class': 'btn-check',
            'autocomplete': 'off'
        }),
        required=True
    )
    xthreshold = forms.IntegerField(min_value=1, 
        required=True,
        initial=20,
        widget=forms.NumberInput(attrs={'class': 'form-control form-control-lg ms-2','style': 'width: 100px;height:40px'}))
    fasta_file = forms.FileField(
        widget=forms.ClearableFileInput(attrs={'class':'form-control','style': 'width: 350px;'}),
        required=True)

    def clean_xthreshold(self):
        xthreshold = self.cleaned_data.get('xthreshold')
        if xthreshold < 1:
            raise forms.ValidationError("O valor mínimo para X-Threshold é 1.")
        return xthreshold