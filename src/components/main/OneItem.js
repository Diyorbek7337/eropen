import React, { useState } from 'react';
import { useTranslation } from 'react-i18next';

const OneItem = () => {
    const { t } = useTranslation()
    const [ block, setBlock ] = useState(false)

    return (
        <div className="block__item">
            <button className={block ? "block__title --active": "block__title"} type="button" data-spoller="" onClick={() => setBlock(!block)}>
            <i className="fas fa-file-medical"></i>{t('Home.compositionTitle')}

            <span className="spoller__arrow"></span>
            </button>
            <div className={block ? "block__text --active": "block__text"}>
            {t('Home.compositionText1')} <br/>
            {t('Home.compositionText2')} <br/>
            {t('Home.compositionText3')}
            </div>
        </div>
    );
};

export default OneItem;